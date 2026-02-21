#include "ikafssnserver/request_processor.hpp"

#include "core/config.hpp"
#include "core/kmer_encoding.hpp"
#include "io/result_writer.hpp"
#include "search/oid_filter.hpp"

#include <algorithm>
#include <mutex>

namespace ikafssn {

SearchResponse process_search_request(
    const SearchRequest& req,
    const std::map<int, KmerGroup>& kmer_groups,
    int default_k,
    const SearchConfig& default_config) {

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
    if (req.min_stage1_score != 0)
        config.stage1.min_stage1_score = req.min_stage1_score;
    if (req.num_results != 0)
        config.num_results = req.num_results;
    if (req.mode != 0)
        config.mode = req.mode;
    if (req.stage1_score_type != 0)
        config.stage1.stage1_score_type = req.stage1_score_type;
    if (req.sort_score != 0)
        config.sort_score = req.sort_score;

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
    resp.mode = config.mode;
    resp.stage1_score_type = config.stage1.stage1_score_type;

    // Build seqidlist filter mode
    OidFilterMode filter_mode = OidFilterMode::kNone;
    if (req.seqidlist_mode == SeqidlistMode::kInclude)
        filter_mode = OidFilterMode::kInclude;
    else if (req.seqidlist_mode == SeqidlistMode::kExclude)
        filter_mode = OidFilterMode::kExclude;

    // Process each query
    resp.status = 0;
    resp.results.resize(req.queries.size());

    for (size_t qi = 0; qi < req.queries.size(); qi++) {
        const auto& query = req.queries[qi];
        auto& qr = resp.results[qi];
        qr.query_id = query.query_id;

        // Collect hits from all volumes
        for (const auto& vol : group.volumes) {
            // Build per-volume OID filter
            OidFilter oid_filter;
            if (filter_mode != OidFilterMode::kNone) {
                oid_filter.build(req.seqids, vol.ksx, filter_mode);
            }

            SearchResult sr;
            if (group.kmer_type == 0) {
                sr = search_volume<uint16_t>(
                    query.query_id, query.sequence, k,
                    vol.kix, vol.kpx, vol.ksx, oid_filter, config);
            } else {
                sr = search_volume<uint32_t>(
                    query.query_id, query.sequence, k,
                    vol.kix, vol.kpx, vol.ksx, oid_filter, config);
            }

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
                qr.hits.push_back(rh);
            }
        }

        if (config.num_results > 0) {
            // Sort by sort_score descending
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

            // Truncate to num_results
            if (qr.hits.size() > config.num_results) {
                qr.hits.resize(config.num_results);
            }
        }
        // num_results == 0: unlimited, skip sort and truncation
    }

    return resp;
}

} // namespace ikafssn
