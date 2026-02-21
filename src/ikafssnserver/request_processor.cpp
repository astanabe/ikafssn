#include "ikafssnserver/request_processor.hpp"
#include "ikafssnserver/server.hpp"

#include "core/config.hpp"
#include "core/kmer_encoding.hpp"
#include <cmath>
#include "io/result_writer.hpp"
#include "search/oid_filter.hpp"
#include "search/query_preprocessor.hpp"

#include <algorithm>

#include <tbb/blocked_range.h>
#include <tbb/combinable.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

namespace ikafssn {

struct AcceptedQuery {
    size_t result_idx;   // index into resp.results
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
    if (req.has_min_score)
        config.stage2.min_score = req.min_score;
    if (req.max_gap != 0)
        config.stage2.max_gap = req.max_gap;
    if (req.chain_max_lookback != 0)
        config.stage2.chain_max_lookback = req.chain_max_lookback;
    if (req.max_freq_frac_x10000 != 0) {
        double frac = static_cast<double>(req.max_freq_frac_x10000) / 10000.0;
        uint64_t total_nseq = 0;
        for (const auto& vol : group.volumes) total_nseq += vol.ksx.num_sequences();
        config.stage1.max_freq = static_cast<uint32_t>(std::ceil(frac * total_nseq));
        if (config.stage1.max_freq == 0) config.stage1.max_freq = 1;
    } else if (req.max_freq != 0) {
        config.stage1.max_freq = req.max_freq;
    }
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
        bool has_min_score = req.has_min_score;
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

    // Resolve effective accept_qdegen: client non-zero overrides server default
    uint8_t accept_qdegen = (req.accept_qdegen != 0)
        ? req.accept_qdegen : config.accept_qdegen;

    // Classify queries: skipped (degenerate) vs valid
    std::vector<size_t> valid_indices;  // original indices of non-skipped queries
    std::vector<size_t> skipped_indices;
    for (size_t qi = 0; qi < req.queries.size(); qi++) {
        if (accept_qdegen == 0 &&
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

    // Build vector of KixReader pointers for global preprocessing
    std::vector<const KixReader*> all_kix;
    all_kix.reserve(group.volumes.size());
    for (const auto& vol : group.volumes) {
        all_kix.push_back(&vol.kix);
    }
    const KhxReader* khx_ptr = group.khx.is_open() ? &group.khx : nullptr;

    // Preprocess accepted queries and build jobs
    struct PreprocessedQuery16 { QueryKmerData<uint16_t> qdata; };
    struct PreprocessedQuery32 { QueryKmerData<uint32_t> qdata; };
    std::vector<PreprocessedQuery16> pp16;
    std::vector<PreprocessedQuery32> pp32;
    // Map from original query index to preprocessed index
    std::vector<size_t> query_pp_idx(req.queries.size(), SIZE_MAX);

    std::vector<AcceptedQuery> accepted_queries;
    accepted_queries.reserve(static_cast<size_t>(acquired));

    // Iterate in original order: include skipped and accepted, skip rejected
    for (size_t qi = 0; qi < req.queries.size(); qi++) {
        bool skipped = (accept_qdegen == 0 &&
                        contains_degenerate_base(req.queries[qi].sequence));
        if (skipped) {
            QueryResult qr;
            qr.query_id = req.queries[qi].query_id;
            qr.skipped = 1;
            resp.results.push_back(std::move(qr));
            continue;
        }

        if (!is_accepted[qi]) continue; // rejected, skip

        // Preprocess this accepted query
        if (group.kmer_type == 0) {
            query_pp_idx[qi] = pp16.size();
            pp16.push_back({preprocess_query<uint16_t>(
                req.queries[qi].sequence, k, all_kix, khx_ptr, config)});
        } else {
            query_pp_idx[qi] = pp32.size();
            pp32.push_back({preprocess_query<uint32_t>(
                req.queries[qi].sequence, k, all_kix, khx_ptr, config)});
        }

        size_t result_idx = resp.results.size();
        QueryResult qr;
        qr.query_id = req.queries[qi].query_id;
        resp.results.push_back(std::move(qr));

        accepted_queries.push_back({result_idx, qi});
    }

    // Thread-local Stage1Buffer to avoid per-job allocation
    uint32_t max_num_seqs = 0;
    for (const auto& vol : group.volumes)
        max_num_seqs = std::max(max_num_seqs, vol.kix.num_sequences());

    tbb::enumerable_thread_specific<Stage1Buffer> tls_bufs(
        [max_num_seqs]() {
            Stage1Buffer buf;
            buf.ensure_capacity(max_num_seqs);
            return buf;
        });

    // Thread-local hit collection: (result_idx, ResponseHit) pairs
    tbb::combinable<std::vector<std::pair<size_t, ResponseHit>>> tls_hits;

    // Parallel execution using preprocessed data (query-level granularity)
    arena.execute([&] {
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, accepted_queries.size()),
            [&](const tbb::blocked_range<size_t>& range) {
                auto& buf = tls_bufs.local();
                auto& local_hits = tls_hits.local();

                for (size_t i = range.begin(); i != range.end(); ++i) {
                    const auto& aq = accepted_queries[i];
                    const auto& query = req.queries[aq.query_idx];
                    size_t pp_idx = query_pp_idx[aq.query_idx];

                    for (size_t vol_i = 0; vol_i < group.volumes.size(); vol_i++) {
                        const auto& vol = group.volumes[vol_i];

                        // Build per-volume OID filter
                        OidFilter oid_filter;
                        if (filter_mode != OidFilterMode::kNone) {
                            oid_filter.build(req.seqids, vol.ksx, filter_mode);
                        }

                        SearchResult sr;
                        if (group.kmer_type == 0) {
                            sr = search_volume<uint16_t>(
                                query.query_id, pp16[pp_idx].qdata, k,
                                vol.kix, vol.kpx, vol.ksx, oid_filter, config, &buf);
                        } else {
                            sr = search_volume<uint32_t>(
                                query.query_id, pp32[pp_idx].qdata, k,
                                vol.kix, vol.kpx, vol.ksx, oid_filter, config, &buf);
                        }

                        if (!sr.hits.empty()) {
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
                                local_hits.emplace_back(aq.result_idx, rh);
                            }
                        }
                    }
                }
            });
    });

    // Merge thread-local hits into resp.results
    tls_hits.combine_each([&resp](std::vector<std::pair<size_t, ResponseHit>>& local) {
        for (auto& [idx, rh] : local) {
            resp.results[idx].hits.push_back(std::move(rh));
        }
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
