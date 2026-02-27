#include "ikafssnserver/request_processor.hpp"
#include "ikafssnserver/server.hpp"

#include "core/config.hpp"
#include "core/kmer_encoding.hpp"
#include <climits>
#include <cmath>
#include "io/fasta_reader.hpp"
#include "io/result_writer.hpp"
#include "search/oid_filter.hpp"
#include "search/query_preprocessor.hpp"

#include <algorithm>
#include <unordered_map>

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
    const DatabaseEntry& db,
    Server& server,
    tbb::task_arena& arena) {

    SearchResponse resp;
    resp.db = db.name;

    // Determine k
    int k = (req.k != 0) ? req.k : db.default_k;
    resp.k = static_cast<uint8_t>(k);

    auto it = db.kmer_groups.find(k);
    if (it == db.kmer_groups.end()) {
        resp.status = 1;
        return resp;
    }

    const KmerGroup& group = it->second;

    // Build search config, using request values if non-zero, else DB defaults
    SearchConfig config = db.resolved_search_config;
    if (req.has_stage2_min_score)
        config.stage2.min_score = req.stage2_min_score;
    if (req.stage2_max_gap != 0)
        config.stage2.max_gap = req.stage2_max_gap;
    if (req.stage2_max_lookback != 0)
        config.stage2.chain_max_lookback = req.stage2_max_lookback;
    if (req.stage1_max_freq_frac_x10000 != 0) {
        double frac = static_cast<double>(req.stage1_max_freq_frac_x10000) / 10000.0;
        uint64_t total_nseq = 0;
        for (const auto& vol : group.volumes) total_nseq += vol.ksx.num_sequences();
        config.stage1.max_freq = static_cast<uint32_t>(std::ceil(frac * total_nseq));
        if (config.stage1.max_freq == 0) config.stage1.max_freq = 1;
    } else if (req.stage1_max_freq != 0) {
        config.stage1.max_freq = req.stage1_max_freq;
    }
    if (req.stage2_min_diag_hits != 0)
        config.stage2.min_diag_hits = req.stage2_min_diag_hits;
    if (req.stage1_topn != 0)
        config.stage1.stage1_topn = req.stage1_topn;
    if (req.stage1_min_score_frac_x10000 != 0) {
        config.min_stage1_score_frac =
            static_cast<double>(req.stage1_min_score_frac_x10000) / 10000.0;
    } else if (req.stage1_min_score != 0) {
        config.stage1.min_stage1_score = req.stage1_min_score;
    }
    if (req.num_results != 0)
        config.num_results = req.num_results;
    if (req.mode != 0)
        config.mode = req.mode;
    if (req.stage1_score != 0)
        config.stage1.stage1_score_type = req.stage1_score;
    if (req.strand != 0)
        config.strand = req.strand;

    // Validate mode against max_mode
    if (config.mode > db.max_mode) {
        resp.status = 4; // mode exceeds max_mode for this DB
        return resp;
    }

    // sort_score is auto-determined by mode
    switch (config.mode) {
        case 1: config.sort_score = 1; break;
        case 2: config.sort_score = 2; break;
        case 3: config.sort_score = 3; break;
        default: config.sort_score = 2; break;
    }

    // Mode 1: consistency checks
    if (config.mode == 1) {
        bool has_min_score = req.has_stage2_min_score;
        bool has_min_stage1_score = (req.stage1_min_score != 0);
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

    // Resolve Stage 3 parameters from request (INT16_MIN = use server default)
    Stage3Config stage3_config = db.stage3_config;
    if (req.stage3_traceback != 0)
        stage3_config.traceback = true;
    if (req.stage3_gapopen != INT16_MIN)
        stage3_config.gapopen = req.stage3_gapopen;
    if (req.stage3_gapext != INT16_MIN)
        stage3_config.gapext = req.stage3_gapext;
    if (req.stage3_min_pident_x100 != 0)
        stage3_config.min_pident = static_cast<double>(req.stage3_min_pident_x100) / 100.0;
    if (req.stage3_min_nident != 0)
        stage3_config.min_nident = req.stage3_min_nident;

    // Resolve context
    bool ctx_is_ratio = db.context_is_ratio;
    double ctx_ratio = db.context_ratio;
    uint32_t ctx_abs = db.context_abs;
    if (req.context_frac_x10000 != 0) {
        ctx_is_ratio = true;
        ctx_ratio = static_cast<double>(req.context_frac_x10000) / 10000.0;
    } else if (req.context_abs != 0) {
        ctx_is_ratio = false;
        ctx_abs = req.context_abs;
    }

    // Set response metadata
    resp.status = 0;
    resp.mode = config.mode;
    resp.stage1_score = config.stage1.stage1_score_type;

    // Build seqidlist filter mode
    OidFilterMode filter_mode = OidFilterMode::kNone;
    if (req.seqidlist_mode == SeqidlistMode::kInclude)
        filter_mode = OidFilterMode::kInclude;
    else if (req.seqidlist_mode == SeqidlistMode::kExclude)
        filter_mode = OidFilterMode::kExclude;

    uint8_t accept_qdegen = req.accept_qdegen;

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
        resp.rejected_qseqids.push_back(req.queries[valid_indices[i]].qseqid);
    }

    // Build resp.results for skipped + accepted queries only
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
    std::vector<size_t> query_pp_idx(req.queries.size(), SIZE_MAX);

    std::vector<AcceptedQuery> accepted_queries;
    accepted_queries.reserve(static_cast<size_t>(acquired));

    // Iterate in original order: include skipped and accepted, skip rejected
    for (size_t qi = 0; qi < req.queries.size(); qi++) {
        bool skipped = (accept_qdegen == 0 &&
                        contains_degenerate_base(req.queries[qi].sequence));
        if (skipped) {
            QueryResult qr;
            qr.qseqid = req.queries[qi].qseqid;
            qr.skipped = 1;
            resp.results.push_back(std::move(qr));
            continue;
        }

        if (!is_accepted[qi]) continue; // rejected, skip

        // Preprocess this accepted query
        bool multi_degen = false;
        if (group.kmer_type == 0) {
            query_pp_idx[qi] = pp16.size();
            pp16.push_back({preprocess_query<uint16_t>(
                req.queries[qi].sequence, k, all_kix, khx_ptr, config)});
            multi_degen = pp16.back().qdata.has_multi_degen;
        } else {
            query_pp_idx[qi] = pp32.size();
            pp32.push_back({preprocess_query<uint32_t>(
                req.queries[qi].sequence, k, all_kix, khx_ptr, config)});
            multi_degen = pp32.back().qdata.has_multi_degen;
        }

        size_t result_idx = resp.results.size();
        QueryResult qr;
        qr.qseqid = req.queries[qi].qseqid;
        if (multi_degen) {
            qr.warnings |= kWarnMultiDegen;
        }
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
                                query.qseqid, pp16[pp_idx].qdata, k,
                                vol.kix, vol.kpx, vol.ksx, oid_filter, config, &buf);
                        } else {
                            sr = search_volume<uint32_t>(
                                query.qseqid, pp32[pp_idx].qdata, k,
                                vol.kix, vol.kpx, vol.ksx, oid_filter, config, &buf);
                        }

                        if (!sr.hits.empty()) {
                            for (const auto& cr : sr.hits) {
                                ResponseHit rh;
                                rh.sseqid = std::string(vol.ksx.accession(cr.seq_id));
                                rh.sstrand = cr.is_reverse ? 1 : 0;
                                rh.qstart = cr.q_start;
                                rh.qend = cr.q_end;
                                rh.sstart = cr.s_start;
                                rh.send = cr.s_end;
                                rh.chainscore = static_cast<uint16_t>(cr.chainscore);
                                if (config.stage1.stage1_score_type == 2)
                                    rh.matchscore = static_cast<uint16_t>(cr.stage1_score);
                                else
                                    rh.coverscore = static_cast<uint16_t>(cr.stage1_score);
                                rh.volume = vol.volume_index;
                                rh.qlen = static_cast<uint32_t>(query.sequence.size());
                                rh.slen = vol.ksx.seq_length(cr.seq_id);
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

    // Stage 3 alignment (mode 3 only)
    if (config.mode == 3) {
        if (db.db_path.empty()) {
            resp.status = 3; // mode 3 requires BLAST DB
            server.release_sequences(acquired);
            return resp;
        }

        // Build FastaRecord vector from accepted queries
        std::vector<FastaRecord> fasta_queries;
        for (const auto& aq : accepted_queries) {
            fasta_queries.push_back({
                req.queries[aq.query_idx].qseqid,
                req.queries[aq.query_idx].sequence
            });
        }

        // Flatten all hits into OutputHit vector for run_stage3
        std::vector<OutputHit> output_hits;
        for (auto& qr : resp.results) {
            for (auto& hit : qr.hits) {
                OutputHit oh;
                oh.qseqid = qr.qseqid;
                oh.sseqid = hit.sseqid;
                oh.sstrand = (hit.sstrand == 0) ? '+' : '-';
                oh.qstart = hit.qstart;
                oh.qend = hit.qend;
                oh.sstart = hit.sstart;
                oh.send = hit.send;
                oh.chainscore = hit.chainscore;
                oh.coverscore = hit.coverscore;
                oh.matchscore = hit.matchscore;
                oh.volume = hit.volume;
                oh.qlen = hit.qlen;
                oh.slen = hit.slen;
                output_hits.push_back(std::move(oh));
            }
        }

        Logger logger(Logger::kInfo);
        output_hits = run_stage3(output_hits, fasta_queries, db.db_path,
                                 stage3_config, ctx_is_ratio, ctx_ratio, ctx_abs, logger);

        // Write back to ResponseHit
        for (auto& qr : resp.results) qr.hits.clear();
        std::unordered_map<std::string, size_t> qid_to_ridx;
        for (size_t i = 0; i < resp.results.size(); i++)
            qid_to_ridx[resp.results[i].qseqid] = i;
        for (const auto& oh : output_hits) {
            auto qit = qid_to_ridx.find(oh.qseqid);
            if (qit == qid_to_ridx.end()) continue;
            ResponseHit rh;
            rh.sseqid = oh.sseqid;
            rh.sstrand = (oh.sstrand == '+') ? 0 : 1;
            rh.qstart = oh.qstart;
            rh.qend = oh.qend;
            rh.sstart = oh.sstart;
            rh.send = oh.send;
            rh.chainscore = static_cast<uint16_t>(oh.chainscore);
            rh.coverscore = static_cast<uint16_t>(oh.coverscore);
            rh.matchscore = static_cast<uint16_t>(oh.matchscore);
            rh.volume = oh.volume;
            rh.qlen = oh.qlen;
            rh.slen = oh.slen;
            rh.alnscore = oh.alnscore;
            rh.nident = oh.nident;
            rh.mismatch = oh.mismatch;
            rh.pident_x100 = static_cast<uint16_t>(oh.pident * 100.0);
            rh.cigar = oh.cigar;
            rh.qseq = oh.qseq;
            rh.sseq = oh.sseq;
            resp.results[qit->second].hits.push_back(std::move(rh));
        }

        resp.stage3_traceback = stage3_config.traceback ? 1 : 0;
    }

    // Release permits
    server.release_sequences(acquired);

    // Post-process: sort/truncate per accepted query
    if (config.num_results > 0) {
        for (auto& qr : resp.results) {
            if (qr.skipped != 0) continue;

            if (config.sort_score == 1) {
                std::sort(qr.hits.begin(), qr.hits.end(),
                          [](const ResponseHit& a, const ResponseHit& b) {
                              return (a.coverscore + a.matchscore) > (b.coverscore + b.matchscore);
                          });
            } else if (config.sort_score == 3) {
                std::sort(qr.hits.begin(), qr.hits.end(),
                          [](const ResponseHit& a, const ResponseHit& b) {
                              return a.alnscore > b.alnscore;
                          });
            } else {
                std::sort(qr.hits.begin(), qr.hits.end(),
                          [](const ResponseHit& a, const ResponseHit& b) {
                              return a.chainscore > b.chainscore;
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
