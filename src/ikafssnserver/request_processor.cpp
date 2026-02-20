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
        std::vector<OutputHit> all_hits;

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
                rh.volume = vol.volume_index;
                qr.hits.push_back(rh);
            }
        }

        // Sort by score descending
        std::sort(qr.hits.begin(), qr.hits.end(),
                  [](const ResponseHit& a, const ResponseHit& b) {
                      return a.score > b.score;
                  });

        // Truncate to num_results
        uint32_t num_results = config.num_results;
        if (qr.hits.size() > num_results) {
            qr.hits.resize(num_results);
        }
    }

    return resp;
}

} // namespace ikafssn
