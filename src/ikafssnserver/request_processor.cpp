#include "ikafssnserver/request_processor.hpp"
#include "ikafssnserver/server.hpp"

#include "core/config.hpp"
#include "core/kmer_encoding.hpp"
#include "io/result_writer.hpp"
#include "search/oid_filter.hpp"

#include <algorithm>
#include <mutex>

#include <tbb/parallel_for_each.h>

namespace ikafssn {

struct SearchJob {
    size_t result_idx;   // index into resp.results
    size_t volume_idx;   // index into group.volumes
    size_t query_idx;    // index into req.queries (for sequence data)
};

SearchResponse process_search_request(
    const SearchRequest& req,
    const std::map<int, KmerGroup>& kmer_groups,
    int default_k,
    const SearchConfig& default_config,
    Server& server,
    tbb::task_arena& arena) {

    SearchResponse resp;

    // Determine k
    int k = (req.k != 0) ? req.k : default_k;
    resp.k = static_cast<uint8_t>(k);

    auto it = kmer_groups.find(k);
    if (it == kmer_groups.end()) {
        resp.status = 1;
        return resp;
    }

    const KmerGroup& group = it->second;

    // Build search config, using request values if non-zero, else server defaults
    SearchConfig config = default_config;
    if (req.min_score != 0)
        config.stage2.min_score = req.min_score;
    if (req.max_gap != 0)
        config.stage2.max_gap = req.max_gap;
    if (req.max_freq != 0)
        config.stage1.max_freq = req.max_freq;
    if (req.min_diag_hits != 0)
        config.stage2.min_diag_hits = req.min_diag_hits;
    if (req.stage1_topn != 0)
        config.stage1.stage1_topn = req.stage1_topn;
    if (req.min_stage1_score_frac_x10000 != 0) {
        config.min_stage1_score_frac =
            static_cast<double>(req.min_stage1_score_frac_x10000) / 10000.0;
    } else if (req.min_stage1_score != 0) {
        config.stage1.min_stage1_score = req.min_stage1_score;
    }
    if (req.num_results != 0)
        config.num_results = req.num_results;
    if (req.mode != 0)
        config.mode = req.mode;
    if (req.stage1_score_type != 0)
        config.stage1.stage1_score_type = req.stage1_score_type;
    if (req.sort_score != 0)
        config.sort_score = req.sort_score;
    if (req.strand != 0)
        config.strand = req.strand;

    // Mode 1: force sort_score=1
    if (config.mode == 1) {
        config.sort_score = 1;

        // Consistency check: min_score and min_stage1_score
        bool has_min_score = (req.min_score != 0);
        bool has_min_stage1_score = (req.min_stage1_score != 0);
        if (has_min_score && has_min_stage1_score &&
            config.stage2.min_score != config.stage1.min_stage1_score) {
            resp.status = 2; // parameter conflict
            return resp;
        }
        if (has_min_score && !has_min_stage1_score) {
            config.stage1.min_stage1_score = config.stage2.min_score;
        }
        if (!has_min_score && has_min_stage1_score) {
            config.stage2.min_score = config.stage1.min_stage1_score;
        }
    }

    // Set response metadata
    resp.status = 0;
    resp.mode = config.mode;
    resp.stage1_score_type = config.stage1.stage1_score_type;

    // Build seqidlist filter mode
    OidFilterMode filter_mode = OidFilterMode::kNone;
    if (req.seqidlist_mode == SeqidlistMode::kInclude)
        filter_mode = OidFilterMode::kInclude;
    else if (req.seqidlist_mode == SeqidlistMode::kExclude)
        filter_mode = OidFilterMode::kExclude;

    // Classify queries: skipped (degenerate) vs valid
    std::vector<size_t> valid_indices;  // original indices of non-skipped queries
    std::vector<size_t> skipped_indices;
    for (size_t qi = 0; qi < req.queries.size(); qi++) {
        if (req.accept_qdegen == 0 &&
            contains_degenerate_base(req.queries[qi].sequence)) {
            skipped_indices.push_back(qi);
        } else {
            valid_indices.push_back(qi);
        }
    }

    int valid_count = static_cast<int>(valid_indices.size());

    // Acquire per-sequence permits
    int acquired = server.try_acquire_sequences(valid_count);

    // First `acquired` valid queries are accepted; rest are rejected
    for (int i = acquired; i < valid_count; i++) {
        resp.rejected_query_ids.push_back(req.queries[valid_indices[i]].query_id);
    }

    // Build resp.results for skipped + accepted queries only
    // (rejected queries are NOT included in results)
    // Mark accepted queries for O(1) lookup
    std::vector<bool> is_accepted(req.queries.size(), false);
    for (int i = 0; i < acquired; i++) {
        is_accepted[valid_indices[i]] = true;
    }

    std::vector<SearchJob> jobs;
    jobs.reserve(static_cast<size_t>(acquired) * group.volumes.size());

    // Iterate in original order: include skipped and accepted, skip rejected
    for (size_t qi = 0; qi < req.queries.size(); qi++) {
        bool skipped = (req.accept_qdegen == 0 &&
                        contains_degenerate_base(req.queries[qi].sequence));
        if (skipped) {
            QueryResult qr;
            qr.query_id = req.queries[qi].query_id;
            qr.skipped = 1;
            resp.results.push_back(std::move(qr));
            continue;
        }

        if (!is_accepted[qi]) continue; // rejected, skip

        size_t result_idx = resp.results.size();
        QueryResult qr;
        qr.query_id = req.queries[qi].query_id;
        resp.results.push_back(std::move(qr));

        // Create jobs for this accepted query
        for (size_t vol_i = 0; vol_i < group.volumes.size(); vol_i++) {
            jobs.push_back({result_idx, vol_i, qi});
        }
    }

    // Parallel execution
    std::mutex hits_mutex;
    arena.execute([&] {
        tbb::parallel_for_each(jobs.begin(), jobs.end(),
            [&](const SearchJob& job) {
                const auto& query = req.queries[job.query_idx];
                const auto& vol = group.volumes[job.volume_idx];

                // Build per-volume OID filter
                OidFilter oid_filter;
                if (filter_mode != OidFilterMode::kNone) {
                    oid_filter.build(req.seqids, vol.ksx, filter_mode);
                }

                const KhxReader* khx_ptr = vol.khx.is_open() ? &vol.khx : nullptr;

                SearchResult sr;
                if (group.kmer_type == 0) {
                    sr = search_volume<uint16_t>(
                        query.query_id, query.sequence, k,
                        vol.kix, vol.kpx, vol.ksx, oid_filter, config, khx_ptr);
                } else {
                    sr = search_volume<uint32_t>(
                        query.query_id, query.sequence, k,
                        vol.kix, vol.kpx, vol.ksx, oid_filter, config, khx_ptr);
                }

                if (!sr.hits.empty()) {
                    std::vector<ResponseHit> local_hits;
                    local_hits.reserve(sr.hits.size());
                    for (const auto& cr : sr.hits) {
                        ResponseHit rh;
                        rh.accession = std::string(vol.ksx.accession(cr.seq_id));
                        rh.strand = cr.is_reverse ? 1 : 0;
                        rh.q_start = cr.q_start;
                        rh.q_end = cr.q_end;
                        rh.s_start = cr.s_start;
                        rh.s_end = cr.s_end;
                        rh.score = static_cast<uint16_t>(cr.score);
                        rh.stage1_score = static_cast<uint16_t>(cr.stage1_score);
                        rh.volume = vol.volume_index;
                        local_hits.push_back(rh);
                    }

                    std::lock_guard<std::mutex> lock(hits_mutex);
                    auto& hits = resp.results[job.result_idx].hits;
                    hits.insert(hits.end(), local_hits.begin(), local_hits.end());
                }
            });
    });

    // Release permits
    server.release_sequences(acquired);

    // Post-process: sort/truncate per accepted query
    if (config.num_results > 0) {
        for (auto& qr : resp.results) {
            if (qr.skipped != 0) continue;

            if (config.sort_score == 1) {
                std::sort(qr.hits.begin(), qr.hits.end(),
                          [](const ResponseHit& a, const ResponseHit& b) {
                              return a.stage1_score > b.stage1_score;
                          });
            } else {
                std::sort(qr.hits.begin(), qr.hits.end(),
                          [](const ResponseHit& a, const ResponseHit& b) {
                              return a.score > b.score;
                          });
            }

            if (qr.hits.size() > config.num_results) {
                qr.hits.resize(config.num_results);
            }
        }
    }

    return resp;
}

} // namespace ikafssn
