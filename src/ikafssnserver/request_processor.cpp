#include "ikafssnserver/request_processor.hpp"
#include "ikafssnserver/server.hpp"

#include "core/config.hpp"
#include "core/kmer_encoding.hpp"
#include "core/spaced_seed.hpp"
#include <climits>
#include <cmath>
#include "io/fasta_reader.hpp"
#include "io/result_writer.hpp"
#include "search/oid_filter.hpp"
#include "search/parallel_search.hpp"
#include "search/preprocess_runner.hpp"
#include "search/query_preprocessor.hpp"
#include "search/tier_selection.hpp"

#include <algorithm>
#include <functional>
#include <type_traits>
#include <unordered_map>
#include <utility>

#include <tbb/blocked_range.h>
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

    // Determine k and t
    int k = (req.k != 0) ? req.k : db.default_k;
    uint8_t t = (req.t != 0) ? req.t : db.default_t;
    uint8_t tt = (req.template_type != 0) ? req.template_type : db.default_template_type;
    resp.k = static_cast<uint8_t>(k);
    resp.t = t;

    // Determine group(s) based on template type
    const KmerGroup* group_ptr = nullptr;
    const KmerGroup* group_cod = nullptr;
    const KmerGroup* group_opt = nullptr;
    bool is_both_mode = (tt == 3);

    if (is_both_mode) {
        auto [cod, opt] = db.find_both_groups(k, t);
        if (!cod || !opt) {
            resp.status = 1;
            return resp;
        }
        if (cod->volumes.size() != opt->volumes.size()) {
            resp.status = 1;
            return resp;
        }
        group_cod = cod;
        group_opt = opt;
        group_ptr = cod; // use coding group for metadata
    } else {
        group_ptr = db.find_group(k, t, tt);
        if (!group_ptr) {
            resp.status = 1;
            return resp;
        }
    }
    const KmerGroup& group = *group_ptr;

    // Build search config, using request values if non-zero, else DB defaults
    SearchConfig config = db.resolved_search_config;
    if (req.has_stage2_min_score)
        config.stage2.min_score = req.stage2_min_score;
    if (req.stage2_max_gap != 0)
        config.stage2.max_gap = req.stage2_max_gap;
    if (req.stage2_max_lookback != 0)
        config.stage2.chain_max_lookback = req.stage2_max_lookback;
    if (req.stage2_max_nhit_per_subject != 0)
        config.stage2.max_nhit_per_subject = req.stage2_max_nhit_per_subject;
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
    if (req.max_degen_expand != 0)
        config.max_degen_expand = req.max_degen_expand;
    config.t = t;
    config.accept_qdegen = req.accept_qdegen;

    // Resolve seed masks for spaced seed preprocessing.
    std::vector<uint32_t> seed_masks;
    std::vector<uint32_t> seed_masks_cod;
    std::vector<uint32_t> seed_masks_opt;
    if (t > 0) {
        if (is_both_mode) {
            seed_masks_cod = get_seed_masks(k, t, TemplateType::kCoding);
            seed_masks_opt = get_seed_masks(k, t, TemplateType::kOptimal);
        } else {
            seed_masks = get_seed_masks(k, t, static_cast<TemplateType>(tt));
        }
    }

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
    if (req.stage3_min_ppositive_x100 != 0)
        stage3_config.min_ppositive = static_cast<double>(req.stage3_min_ppositive_x100) / 100.0;
    if (req.stage3_min_npositive != 0)
        stage3_config.min_npositive = req.stage3_min_npositive;
    if (req.score_matrix != 0) {
        switch (req.score_matrix) {
            case 1: stage3_config.score_matrix = "degmatch"; break;
            case 2: stage3_config.score_matrix = "dnafull"; break;
            case 3: stage3_config.score_matrix = "nuc44"; break;
            default: break;
        }
    }

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

    // All queries are valid for permit acquisition; skip detection happens
    // inside preprocess_query and is reflected in QueryResult.skip_reason.
    int valid_count = static_cast<int>(req.queries.size());

    // Acquire per-sequence permits
    int acquired = server.try_acquire_sequences(valid_count);

    // First `acquired` queries are accepted; rest are rejected
    for (int i = acquired; i < valid_count; i++) {
        resp.rejected_qseqids.push_back(req.queries[i].qseqid);
    }

    // Build resp.results for accepted queries only
    std::vector<bool> is_accepted(req.queries.size(), false);
    for (int i = 0; i < acquired; i++) {
        is_accepted[i] = true;
    }

    // KhxReader pointers (used for fractional Stage 1 threshold's Nhighfreq)
    const KhxReader* khx_ptr = nullptr;
    const KhxReader* khx_ptr_cod = nullptr;
    const KhxReader* khx_ptr_opt = nullptr;

    if (is_both_mode) {
        khx_ptr_cod = group_cod->khx.is_open() ? &group_cod->khx : nullptr;
        khx_ptr_opt = group_opt->khx.is_open() ? &group_opt->khx : nullptr;
    } else {
        khx_ptr = group.khx.is_open() ? &group.khx : nullptr;
    }

    // Preprocess accepted queries in parallel and build jobs.
    // Slots are pre-allocated to req.queries.size(); rejected slots remain
    // default-constructed and are simply not referenced downstream.
    struct PreprocessedQuery16 { QueryKmerData<uint16_t> qdata; };
    struct PreprocessedQuery32 { QueryKmerData<uint32_t> qdata; };
    std::vector<PreprocessedQuery16> pp16;
    std::vector<PreprocessedQuery32> pp32;
    std::vector<PreprocessedQuery16> pp16_cod;
    std::vector<PreprocessedQuery16> pp16_opt;
    std::vector<PreprocessedQuery32> pp32_cod;
    std::vector<PreprocessedQuery32> pp32_opt;

    if (is_both_mode) {
        if (group.kmer_type == 0) {
            pp16_cod.resize(req.queries.size());
            pp16_opt.resize(req.queries.size());
        } else {
            pp32_cod.resize(req.queries.size());
            pp32_opt.resize(req.queries.size());
        }
    } else {
        if (group.kmer_type == 0) pp16.resize(req.queries.size());
        else                      pp32.resize(req.queries.size());
    }

    // qi -> qi identity mapping after pre-allocation.
    std::vector<size_t> query_pp_idx(req.queries.size());
    for (size_t qi = 0; qi < req.queries.size(); ++qi) query_pp_idx[qi] = qi;

    // Per-query skip / multi-degen state filled inside the parallel body.
    std::vector<uint8_t> per_q_skip_reason(req.queries.size(), 0);
    std::vector<std::string> per_q_skip_detail(req.queries.size());
    std::vector<uint8_t> per_q_multi_degen(req.queries.size(), 0);

    arena.execute([&] {
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, req.queries.size()),
            [&](const tbb::blocked_range<size_t>& range) {
                for (size_t qi = range.begin(); qi != range.end(); ++qi) {
                    if (!is_accepted[qi]) continue;
                    if (group.kmer_type == 0) {
                        std::vector<PreprocessTemplateBinding<uint16_t>> b;
                        if (is_both_mode) {
                            b.resize(2);
                            b[0].slot  = &pp16_cod[qi].qdata;
                            b[0].khx   = khx_ptr_cod;
                            b[0].masks = &seed_masks_cod;
                            b[1].slot  = &pp16_opt[qi].qdata;
                            b[1].khx   = khx_ptr_opt;
                            b[1].masks = &seed_masks_opt;
                        } else {
                            b.resize(1);
                            b[0].slot  = &pp16[qi].qdata;
                            b[0].khx   = khx_ptr;
                            b[0].masks = &seed_masks;
                        }
                        auto out = run_preprocess_one_query<uint16_t>(
                            req.queries[qi].sequence, k, t, config,
                            b.data(), b.size());
                        per_q_skip_reason[qi] = out.skip_reason;
                        per_q_skip_detail[qi] = std::move(out.skip_detail);
                        per_q_multi_degen[qi] = out.multi_degen ? 1 : 0;
                    } else {
                        std::vector<PreprocessTemplateBinding<uint32_t>> b;
                        if (is_both_mode) {
                            b.resize(2);
                            b[0].slot  = &pp32_cod[qi].qdata;
                            b[0].khx   = khx_ptr_cod;
                            b[0].masks = &seed_masks_cod;
                            b[1].slot  = &pp32_opt[qi].qdata;
                            b[1].khx   = khx_ptr_opt;
                            b[1].masks = &seed_masks_opt;
                        } else {
                            b.resize(1);
                            b[0].slot  = &pp32[qi].qdata;
                            b[0].khx   = khx_ptr;
                            b[0].masks = &seed_masks;
                        }
                        auto out = run_preprocess_one_query<uint32_t>(
                            req.queries[qi].sequence, k, t, config,
                            b.data(), b.size());
                        per_q_skip_reason[qi] = out.skip_reason;
                        per_q_skip_detail[qi] = std::move(out.skip_detail);
                        per_q_multi_degen[qi] = out.multi_degen ? 1 : 0;
                    }
                }
            });
    });

    // Sequential second pass: build resp.results in original query order so
    // that result_idx and accepted_queries match the legacy layout, and
    // adjust permit bookkeeping for skipped queries.
    std::vector<AcceptedQuery> accepted_queries;
    accepted_queries.reserve(static_cast<size_t>(acquired));
    for (size_t qi = 0; qi < req.queries.size(); ++qi) {
        if (!is_accepted[qi]) continue;
        size_t result_idx = resp.results.size();
        QueryResult qr;
        qr.qseqid = req.queries[qi].qseqid;
        qr.qlen = static_cast<uint32_t>(req.queries[qi].sequence.size());
        qr.skip_reason = per_q_skip_reason[qi];
        qr.skip_detail = std::move(per_q_skip_detail[qi]);
        if (per_q_multi_degen[qi]) {
            qr.warnings |= kWarnMultiDegen;
        }
        resp.results.push_back(std::move(qr));
        if (per_q_skip_reason[qi] == 0) {
            accepted_queries.push_back({result_idx, qi});
        } else {
            // Skipped queries do not run Stage 1/2/3.  Release the permit so
            // queue_depth reflects actual concurrent work and the next
            // request can claim the slot.
            server.release_sequences(1);
            acquired--;
        }
    }

    // Thread-local Stage1Buffer to avoid per-job allocation
    uint32_t max_num_seqs = 0;
    if (is_both_mode) {
        for (const auto& vol : group_cod->volumes)
            max_num_seqs = std::max(max_num_seqs, vol.kix.num_sequences());
    } else {
        for (const auto& vol : group.volumes)
            max_num_seqs = std::max(max_num_seqs, vol.kix.num_sequences());
    }

    // Determine optimal tier from actual preprocessed k-mer counts.
    uint32_t max_kmer_positions = 0;
    if (is_both_mode) {
        max_kmer_positions = accumulate_max_kmer_positions(max_kmer_positions, pp16_cod);
        max_kmer_positions = accumulate_max_kmer_positions(max_kmer_positions, pp16_opt);
        max_kmer_positions = accumulate_max_kmer_positions(max_kmer_positions, pp32_cod);
        max_kmer_positions = accumulate_max_kmer_positions(max_kmer_positions, pp32_opt);
    } else {
        max_kmer_positions = accumulate_max_kmer_positions(max_kmer_positions, pp16);
        max_kmer_positions = accumulate_max_kmer_positions(max_kmer_positions, pp32);
    }
    Stage1Tier tier = select_tier(max_kmer_positions, max_kmer_positions);

    tbb::enumerable_thread_specific<Stage1Buffer> tls_bufs(
        [max_num_seqs, tier]() {
            Stage1Buffer buf;
            buf.tier = tier;
            buf.ensure_capacity(max_num_seqs);
            return buf;
        });

    // Number of volumes for the search loop
    size_t num_volumes = is_both_mode ? group_cod->volumes.size() : group.volumes.size();

    // Hoist OidFilter construction out of the inner (query, volume) loop:
    // req.seqids, ksx, and filter_mode are all request-invariant, so M*V
    // redundant builds collapse to V builds.  Each OidFilter::build() walks
    // the whole ksx accession array, so parallel_for over volumes amortises
    // the cost.
    std::vector<OidFilter> oid_filters(num_volumes);
    if (filter_mode != OidFilterMode::kNone) {
        arena.execute([&] {
            tbb::parallel_for(
                tbb::blocked_range<size_t>(0, num_volumes),
                [&](const tbb::blocked_range<size_t>& range) {
                    for (size_t vi = range.begin(); vi != range.end(); ++vi) {
                        const auto& vol = is_both_mode
                            ? group_cod->volumes[vi]
                            : group.volumes[vi];
                        oid_filters[vi].build(req.seqids, vol.ksx, filter_mode);
                    }
                });
        });
    }

    // Map each accepted query's req-index into its result-index so we can
    // attach orchestrator hits back to the correct QueryResult slot.
    std::vector<size_t> query_to_result(req.queries.size(), SIZE_MAX);
    for (const auto& aq : accepted_queries) {
        query_to_result[aq.query_idx] = aq.result_idx;
    }

    auto run_orchestrated = [&](auto kmer_int_tag) {
        using KmerInt = decltype(kmer_int_tag);

        std::vector<VolumeBundle<KmerInt>> volume_bundles(num_volumes);
        for (size_t vi = 0; vi < num_volumes; ++vi) {
            const auto& vol_primary = is_both_mode
                ? group_cod->volumes[vi]
                : group.volumes[vi];
            volume_bundles[vi].kix    = &vol_primary.kix;
            volume_bundles[vi].kpx    = &vol_primary.kpx;
            volume_bundles[vi].ksx    = &vol_primary.ksx;
            volume_bundles[vi].filter = &oid_filters[vi];
            volume_bundles[vi].volume_index = vol_primary.volume_index;
            if (is_both_mode) {
                volume_bundles[vi].kix_opt = &group_opt->volumes[vi].kix;
                volume_bundles[vi].kpx_opt = &group_opt->volumes[vi].kpx;
            }
        }

        std::vector<QueryBundle<KmerInt>> query_bundles(req.queries.size());
        std::vector<std::pair<size_t, size_t>> jobs;
        jobs.reserve(accepted_queries.size() * num_volumes);
        for (const auto& aq : accepted_queries) {
            size_t qi = aq.query_idx;
            size_t pp_idx = query_pp_idx[qi];
            query_bundles[qi].query_id = &req.queries[qi].qseqid;
            if (is_both_mode) {
                if constexpr (std::is_same_v<KmerInt, uint16_t>) {
                    query_bundles[qi].qdata_primary   = &pp16_cod[pp_idx].qdata;
                    query_bundles[qi].qdata_secondary = &pp16_opt[pp_idx].qdata;
                } else {
                    query_bundles[qi].qdata_primary   = &pp32_cod[pp_idx].qdata;
                    query_bundles[qi].qdata_secondary = &pp32_opt[pp_idx].qdata;
                }
            } else {
                if constexpr (std::is_same_v<KmerInt, uint16_t>) {
                    query_bundles[qi].qdata_primary = &pp16[pp_idx].qdata;
                } else {
                    query_bundles[qi].qdata_primary = &pp32[pp_idx].qdata;
                }
            }
            for (size_t vi = 0; vi < num_volumes; ++vi) {
                jobs.emplace_back(qi, vi);
            }
        }

        return run_search_jobs<KmerInt>(query_bundles, volume_bundles, jobs,
                                        k, config, is_both_mode, arena, tls_bufs);
    };

    std::vector<OrchestratorHit> orch_hits;
    if (group.kmer_type == 0) orch_hits = run_orchestrated(uint16_t{});
    else                      orch_hits = run_orchestrated(uint32_t{});

    // Convert OrchestratorHit -> ResponseHit at the boundary, fanning out
    // into the corresponding QueryResult slot.
    for (const auto& oh_in : orch_hits) {
        const auto& cr = oh_in.cr;
        size_t result_idx = query_to_result[oh_in.query_idx];
        if (result_idx == SIZE_MAX) continue;  // defensive; should not happen
        const auto& vol = is_both_mode
            ? group_cod->volumes[oh_in.volume_idx]
            : group.volumes[oh_in.volume_idx];
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
        rh.oid = cr.seq_id;
        rh.qlen = static_cast<uint32_t>(req.queries[oh_in.query_idx].sequence.size());
        rh.slen = vol.ksx.seq_length(cr.seq_id);
        resp.results[result_idx].hits.push_back(std::move(rh));
    }

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
                oh.oid = hit.oid;
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
            rh.npositive = oh.npositive;
            rh.nnegative = oh.nnegative;
            rh.ppositive_x100 = static_cast<uint16_t>(oh.ppositive * 100.0);
            rh.cigar = oh.cigar;
            rh.qseq = oh.qseq;
            rh.sseq = oh.sseq;
            resp.results[qit->second].hits.push_back(std::move(rh));
        }

        resp.stage3_traceback = stage3_config.traceback ? 1 : 0;
    }

    // Release permits
    server.release_sequences(acquired);

    // Post-process: sort/truncate per accepted query (parallel across queries)
    if (config.num_results > 0) {
        std::function<bool(const ResponseHit&, const ResponseHit&)> cmp;
        if (config.sort_score == 1) {
            cmp = [](const ResponseHit& a, const ResponseHit& b) {
                return (a.coverscore + a.matchscore) > (b.coverscore + b.matchscore);
            };
        } else if (config.sort_score == 3) {
            cmp = [](const ResponseHit& a, const ResponseHit& b) {
                return a.alnscore > b.alnscore;
            };
        } else {
            cmp = [](const ResponseHit& a, const ResponseHit& b) {
                return a.chainscore > b.chainscore;
            };
        }

        tbb::parallel_for(size_t(0), resp.results.size(), [&](size_t ri) {
            auto& qr = resp.results[ri];
            if (qr.skip_reason != 0) return;
            if (qr.hits.size() <= config.num_results) {
                std::sort(qr.hits.begin(), qr.hits.end(), cmp);
                return;
            }
            std::nth_element(qr.hits.begin(),
                             qr.hits.begin() + config.num_results,
                             qr.hits.end(), cmp);
            qr.hits.resize(config.num_results);
            std::sort(qr.hits.begin(), qr.hits.end(), cmp);
        });
    }

    return resp;
}

} // namespace ikafssn
