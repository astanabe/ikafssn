#include "search/search_orchestrator.hpp"

#include "index/kix_reader.hpp"
#include "index/kpx_reader.hpp"
#include "index/ksx_reader.hpp"
#include "search/parallel_search.hpp"
#include "search/oid_filter.hpp"
#include "search/query_preprocessor.hpp"
#include "search/result_dedup.hpp"
#include "search/stage1_filter.hpp"
#include "util/common_init.hpp"
#include "util/logger.hpp"

#include <algorithm>
#include <cstdio>
#include <functional>
#include <unordered_map>
#include <utility>
#include <vector>

#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_scan.h>
#include <tbb/task_arena.h>

namespace ikafssn {

namespace {

// Pin the dictionary (WILLNEED + HUGEPAGE) and leave the posting file
// MADV_RANDOM. Posting access is random and the working set far exceeds RAM at
// scale, so prefetching posting lists only triggers the kernel readahead window
// per access and roughly doubles device reads with no benefit.
void apply_stage1_madvise(KixReader& kix) {
    kix.apply_madvise_dict_only();
}

// Stage 2A (mode 2/3): same rationale, applied to both kix and kpx.
void apply_stage2a_madvise(KixReader& kix, KpxReader& kpx) {
    kix.apply_madvise_dict_only();
    kpx.apply_madvise_dict_only();
}

template <typename KmerInt>
std::vector<ExtJob> build_ext_jobs(const RunSearchInputs<KmerInt>& in,
                                    size_t num_volumes) {
    const auto& queries = *in.queries;
    const auto* skip_reason = in.query_skip_reason;
    const int8_t strand = in.config.strand;

    std::vector<ExtJob> ext_jobs;
    ext_jobs.reserve(queries.size() * num_volumes * (strand == 2 ? 2 : 1));

    for (size_t qi = 0; qi < queries.size(); ++qi) {
        if (skip_reason && qi < skip_reason->size() && (*skip_reason)[qi] != 0)
            continue;
        for (size_t vi = 0; vi < num_volumes; ++vi) {
            if (strand == 2) {
                ext_jobs.push_back({qi, vi, 0});
                ext_jobs.push_back({qi, vi, 1});
            } else if (strand == 1) {
                ext_jobs.push_back({qi, vi, 0});
            } else {  // -1
                ext_jobs.push_back({qi, vi, 1});
            }
        }
    }
    return ext_jobs;
}

// Build the per-batch index list (indices into ext_jobs) for the volumes in
// the batch.  Linear scan is fine because num_volumes is small relative to
// ext_jobs.size().
std::vector<size_t> indices_for_batch(const std::vector<ExtJob>& ext_jobs,
                                       const std::vector<uint16_t>& vols) {
    std::vector<bool> in_batch;
    uint16_t max_v = 0;
    for (uint16_t v : vols) max_v = std::max(max_v, v);
    in_batch.assign(static_cast<size_t>(max_v) + 1, false);
    for (uint16_t v : vols) in_batch[v] = true;

    std::vector<size_t> out;
    out.reserve(ext_jobs.size());
    for (size_t i = 0; i < ext_jobs.size(); ++i) {
        size_t vi = ext_jobs[i].vi;
        if (vi < in_batch.size() && in_batch[vi]) out.push_back(i);
    }
    return out;
}

// L-th highest value in `v` (0 if fewer than L); reorders v.
uint32_t lth_highest(std::vector<uint32_t>& v, uint32_t L) {
    if (L == 0 || v.size() < L) return 0;
    std::nth_element(v.begin(), v.begin() + (L - 1), v.end(),
                     std::greater<uint32_t>());
    return v[L - 1];
}

// In-total (L) limit for mode 1: keep, per query, the top-L mode1 results plus
// every result tying the L-th Stage 1 score (coverscore == cr.stage1_score).
void apply_in_total_mode1(std::vector<OrchestratorHit>& hits, uint32_t L) {
    if (L == 0 || hits.empty()) return;
    std::unordered_map<size_t, std::vector<uint32_t>> by_q;
    for (const auto& h : hits) by_q[h.query_idx].push_back(h.cr.stage1_score);
    std::unordered_map<size_t, uint32_t> thr;
    for (auto& kv : by_q) thr[kv.first] = lth_highest(kv.second, L);

    std::vector<OrchestratorHit> out;
    out.reserve(hits.size());
    for (auto& h : hits) {
        uint32_t t = thr[h.query_idx];
        if (t == 0 || h.cr.stage1_score >= t) out.push_back(std::move(h));
    }
    hits = std::move(out);
}

// In-total (L) limit for Stage 2: keep, per query (across all volumes and
// strands), the top-L chains by chainscore plus every chain tying the L-th
// chainscore.  Applied after the per-subject selector, before Stage 3 consumes
// the chains, so it also bounds the alignment work in mode 3.
void apply_in_total_stage2(std::vector<OrchestratorHit>& hits, uint32_t L) {
    if (L == 0 || hits.empty()) return;
    std::unordered_map<size_t, std::vector<uint32_t>> by_q;
    for (const auto& h : hits) by_q[h.query_idx].push_back(h.cr.chainscore);
    std::unordered_map<size_t, uint32_t> thr;
    for (auto& kv : by_q) thr[kv.first] = lth_highest(kv.second, L);

    std::vector<OrchestratorHit> out;
    out.reserve(hits.size());
    for (auto& h : hits) {
        uint32_t t = thr[h.query_idx];
        if (t == 0 || h.cr.chainscore >= t) out.push_back(std::move(h));
    }
    hits = std::move(out);
}

// In-total (L) limit for modes 2/3: keep, per query (across all volumes and
// strands), the top-L Stage 1 candidates plus every candidate tying the L-th
// coverscore, by pruning each ext_job's candidate list.  Derived per-ext_job
// state is rebuilt to stay consistent with the pruned candidates.
void apply_in_total_stage1(std::vector<JobState>& states,
                           const std::vector<ExtJob>& ext_jobs, uint32_t L) {
    if (L == 0) return;
    size_t max_qi = 0;
    for (const auto& ej : ext_jobs) max_qi = std::max(max_qi, ej.qi);
    std::vector<std::vector<uint32_t>> scores_by_q(max_qi + 1);
    for (size_t ji = 0; ji < ext_jobs.size(); ++ji) {
        for (const auto& c : states[ji].candidates)
            scores_by_q[ext_jobs[ji].qi].push_back(c.score);
    }
    std::vector<uint32_t> thr(max_qi + 1, 0);
    for (size_t qi = 0; qi <= max_qi; ++qi)
        thr[qi] = lth_highest(scores_by_q[qi], L);

    for (size_t ji = 0; ji < ext_jobs.size(); ++ji) {
        uint32_t t = thr[ext_jobs[ji].qi];
        if (t == 0) continue;
        JobState& st = states[ji];
        if (st.candidates.empty()) continue;
        size_t before = st.candidates.size();
        std::vector<Stage1Candidate> kept;
        kept.reserve(st.candidates.size());
        for (const auto& c : st.candidates)
            if (c.score >= t) kept.push_back(c);
        if (kept.size() == before) continue;
        st.candidates = std::move(kept);
        // Rebuild derived structures consumed by Stage 2.
        st.stage1_scores.clear();
        st.stage1_scores.reserve(st.candidates.size());
        for (const auto& c : st.candidates) st.stage1_scores[c.id] = c.score;
        if (st.both_mode) {
            st.sorted_candidate_sids.clear();
            st.sorted_candidate_sids.reserve(st.candidates.size());
            for (const auto& c : st.candidates)
                st.sorted_candidate_sids.push_back(c.id);
            std::sort(st.sorted_candidate_sids.begin(),
                      st.sorted_candidate_sids.end());
        }
    }
}

}  // namespace

template <typename KmerInt>
std::vector<OrchestratorHit> run_search(const RunSearchInputs<KmerInt>& in) {
    const auto& queries = *in.queries;
    const size_t num_volumes = in.volumes_cod.size();

    Logger* logger = in.logger;

    // Build ext_jobs.
    std::vector<ExtJob> ext_jobs = build_ext_jobs<KmerInt>(in, num_volumes);
    if (ext_jobs.empty()) return {};

    std::vector<JobState> states(ext_jobs.size());

    // Cache wire-level volume_index per volume bundle slot (volume_idx).
    std::vector<uint16_t> volume_indices(num_volumes);
    for (size_t vi = 0; vi < num_volumes; ++vi) {
        volume_indices[vi] = in.volumes_cod[vi].volume_index;
    }

    // Per-volume readers; opened/closed per batch.
    std::vector<KixReader> kix_cod(num_volumes), kix_opt;
    std::vector<KpxReader> kpx_cod(num_volumes), kpx_opt;
    if (in.both_mode) {
        kix_opt.resize(num_volumes);
        kpx_opt.resize(num_volumes);
    }

    // Per-volume bundle storage; pointers refilled per batch, identical
    // shape across stages.
    std::vector<VolumeBundle<KmerInt>> bundles(num_volumes);
    for (size_t vi = 0; vi < num_volumes; ++vi) {
        bundles[vi].ksx          = in.ksx_per_volume[vi];
        bundles[vi].filter       = &in.oid_filters[vi];
        bundles[vi].volume_index = in.volumes_cod[vi].volume_index;
    }

    tbb::task_arena arena(in.nthread);
    auto make_tls_buf = [&in]() {
        Stage1Buffer buf;
        buf.width = in.width;
        buf.ensure_capacity(in.max_num_seqs);
        return buf;
    };
    tbb::enumerable_thread_specific<Stage1Buffer> tls_bufs(make_tls_buf);

    // Each volume contributes the same number of ext_jobs (every query is
    // searched against every volume), so a running counter that bundles
    // volumes until their ext_job total reaches nthread saturates the
    // arena: many-query runs get one volume per group (fastest), few-query
    // runs bundle volumes to fill the thread pool.
    const size_t ext_jobs_per_volume =
        num_volumes > 0 ? ext_jobs.size() / num_volumes : 0;
    const size_t group_ext_job_target = static_cast<size_t>(
        in.nthread > 0 ? in.nthread : 1);

    // Mode 1 fast path uses this batch-local accumulator: each Stage 1
    // batch flushes its mode1_results into it and releases the per-ji
    // capacity, so the orchestrator's peak heap stays bounded by one
    // batch's worth of candidates instead of all batches combined.
    std::vector<OrchestratorHit> mode1_results;

    // ----------------------------------------------------------------
    // Stage 1: incremental open + madvise + accumulate.  Volumes are
    // bundled into a group until the group's ext_job count reaches the
    // thread target, then the group is run, folded, and closed.
    // ----------------------------------------------------------------
    std::vector<uint16_t> s1_cur_vols;
    size_t                s1_cur_ext_jobs = 0;
    size_t                s1_batch_seq = 0;

    auto run_stage1_batch_and_close = [&]() {
        if (s1_cur_vols.empty()) return;
        auto idxs = indices_for_batch(ext_jobs, s1_cur_vols);
        run_stage1_jobs<KmerInt>(idxs, ext_jobs, queries, bundles,
                                  in.k, in.config, in.both_mode,
                                  states, arena, tls_bufs);

        // Mode 1: drain this batch's mode1_results into the run
        // accumulator using a 2-pass parallel_scan + parallel_for so
        // the fold scales with the batch and per-ji capacity is freed
        // immediately.
        if (in.config.mode == 1 && !idxs.empty()) {
            const size_t n_batch = idxs.size();
            std::vector<size_t> offsets(n_batch + 1, 0);
            arena.execute([&] {
                tbb::parallel_scan(
                    tbb::blocked_range<size_t>(0, n_batch),
                    size_t{0},
                    [&](const tbb::blocked_range<size_t>& r, size_t acc, bool is_final) {
                        size_t sum = acc;
                        for (size_t i = r.begin(); i != r.end(); ++i) {
                            if (is_final) offsets[i] = sum;
                            sum += states[idxs[i]].mode1_results.size();
                        }
                        if (is_final && r.end() == n_batch) offsets[n_batch] = sum;
                        return sum;
                    },
                    [](size_t a, size_t b) { return a + b; });
            });

            const size_t batch_total = offsets[n_batch];
            if (batch_total > 0) {
                const size_t base = mode1_results.size();
                mode1_results.resize(base + batch_total);
                arena.execute([&] {
                    tbb::parallel_for(
                        tbb::blocked_range<size_t>(0, n_batch),
                        [&](const tbb::blocked_range<size_t>& r) {
                            for (size_t i = r.begin(); i != r.end(); ++i) {
                                const size_t ji = idxs[i];
                                const ExtJob& ej = ext_jobs[ji];
                                JobState& st = states[ji];
                                size_t dst = base + offsets[i];
                                for (auto& cr : st.mode1_results) {
                                    OrchestratorHit& oh = mode1_results[dst++];
                                    oh.query_idx    = ej.qi;
                                    oh.volume_idx   = ej.vi;
                                    oh.volume_index = volume_indices[ej.vi];
                                    oh.cr           = std::move(cr);
                                }
                                std::vector<ChainResult>().swap(st.mode1_results);
                            }
                        });
                });
            } else {
                arena.execute([&] {
                    tbb::parallel_for(
                        tbb::blocked_range<size_t>(0, n_batch),
                        [&](const tbb::blocked_range<size_t>& r) {
                            for (size_t i = r.begin(); i != r.end(); ++i) {
                                std::vector<ChainResult>().swap(
                                    states[idxs[i]].mode1_results);
                            }
                        });
                });
            }
        }

        if (logger) {
            logger->info("Stage 1 batch %zu: %zu vol(s), %zu ext_job(s)",
                         s1_batch_seq, s1_cur_vols.size(), idxs.size());
        }
        ++s1_batch_seq;

        for (uint16_t vi : s1_cur_vols) {
            kix_cod[vi].close();
            bundles[vi].kix = nullptr;
            if (in.both_mode) {
                kix_opt[vi].close();
                bundles[vi].kix_opt = nullptr;
            }
        }
        s1_cur_vols.clear();
        s1_cur_ext_jobs = 0;
    };

    for (size_t vi = 0; vi < num_volumes; ++vi) {
        if (!kix_cod[vi].open(in.volumes_cod[vi].files.kix_path)) {
            if (logger) logger->error("Cannot open %s",
                in.volumes_cod[vi].files.kix_path.c_str());
            return {};
        }
        if (in.both_mode) {
            if (!kix_opt[vi].open(in.volumes_opt[vi].files.kix_path)) {
                if (logger) logger->error("Cannot open %s",
                    in.volumes_opt[vi].files.kix_path.c_str());
                kix_cod[vi].close();
                return {};
            }
        }

        bundles[vi].kix = &kix_cod[vi];
        bundles[vi].kpx = nullptr;
        if (in.both_mode) {
            bundles[vi].kix_opt = &kix_opt[vi];
            bundles[vi].kpx_opt = nullptr;
        }

        apply_stage1_madvise(kix_cod[vi]);
        if (in.both_mode) {
            apply_stage1_madvise(kix_opt[vi]);
        }
        s1_cur_vols.push_back(static_cast<uint16_t>(vi));
        s1_cur_ext_jobs += ext_jobs_per_volume;
        if (s1_cur_ext_jobs >= group_ext_job_target) {
            run_stage1_batch_and_close();
        }
    }
    run_stage1_batch_and_close();

    if (logger) {
        size_t total_candidates = 0;
        for (auto& st : states) total_candidates += st.candidates.size();
        logger->info("Stage 1 complete: %zu candidate(s) across %zu ext_job(s)",
                     total_candidates, ext_jobs.size());
    }

    // ----------------------------------------------------------------
    // Mode 1 fast path: all per-batch folds were drained above.
    // ----------------------------------------------------------------
    if (in.config.mode == 1) {
        apply_in_total_mode1(mode1_results, in.config.stage1.max_nhit_in_total);
        return mode1_results;
    }

    // Stage 1 in-total (L) limit: prune each query's candidate set across all
    // volumes / strands before Stage 2 consumes it.
    apply_in_total_stage1(states, ext_jobs, in.config.stage1.max_nhit_in_total);

    // ----------------------------------------------------------------
    // Stage 2 (Stage 2A + Stage 2B fused into one incremental loop)
    //
    // Stage 2A and Stage 2B share this loop: Stage 2B runs immediately
    // after Stage 2A for the same batch and the per-ext_job transient
    // state is freed before the next batch starts, so peak memory stays
    // bounded by one group rather than the total volume × query fan-out.
    // ----------------------------------------------------------------
    std::vector<OrchestratorHit> results;
    size_t total_stage2a_hits = 0;
    std::vector<uint16_t> s2_cur_vols;
    size_t                s2_cur_ext_jobs = 0;
    size_t                s2_batch_seq = 0;

    auto run_stage2_batch_and_close = [&]() {
        if (s2_cur_vols.empty()) return;
        auto idxs = indices_for_batch(ext_jobs, s2_cur_vols);
        run_stage2a_jobs<KmerInt>(idxs, ext_jobs, queries, bundles,
                                   in.both_mode, states, arena);

        // Close .kix/.kpx readers before Stage 2B — Stage 2B reads only
        // JobState and does not need the posting files.
        for (uint16_t vi : s2_cur_vols) {
            kpx_cod[vi].close();
            kix_cod[vi].close();
            bundles[vi].kix = nullptr;
            bundles[vi].kpx = nullptr;
            if (in.both_mode) {
                kpx_opt[vi].close();
                kix_opt[vi].close();
                bundles[vi].kix_opt = nullptr;
                bundles[vi].kpx_opt = nullptr;
            }
        }

        if (logger) {
            size_t batch_hits = 0;
            for (size_t ji : idxs) {
                for (auto& kv : states[ji].hits_per_seq)
                    batch_hits += kv.second.size();
            }
            total_stage2a_hits += batch_hits;
            logger->info("Stage 2A batch %zu: %zu vol(s), %zu ext_job(s), %zu hit(s)",
                         s2_batch_seq, s2_cur_vols.size(), idxs.size(),
                         batch_hits);
        }

        size_t chains_before = results.size();
        run_stage2b_jobs(idxs, ext_jobs, states, volume_indices, arena, results);
        if (logger) {
            logger->info("Stage 2B batch %zu: %zu chain(s)",
                         s2_batch_seq, results.size() - chains_before);
        }
        ++s2_batch_seq;

        for (size_t ji : idxs) {
            JobState& st = states[ji];
            st.hits_per_seq  = {};
            st.stage1_scores = {};
            std::vector<SeqId>().swap(st.sorted_candidate_sids);
            std::vector<Stage1Candidate>().swap(st.candidates);
        }

        s2_cur_vols.clear();
        s2_cur_ext_jobs = 0;
    };

    for (size_t vi = 0; vi < num_volumes; ++vi) {
        if (!kix_cod[vi].open(in.volumes_cod[vi].files.kix_path)) {
            if (logger) logger->error("Cannot open %s",
                in.volumes_cod[vi].files.kix_path.c_str());
            return {};
        }
        if (!kpx_cod[vi].open(in.volumes_cod[vi].files.kpx_path)) {
            if (logger) logger->error("Cannot open %s",
                in.volumes_cod[vi].files.kpx_path.c_str());
            kix_cod[vi].close();
            return {};
        }
        if (in.both_mode) {
            if (!kix_opt[vi].open(in.volumes_opt[vi].files.kix_path)) {
                if (logger) logger->error("Cannot open %s",
                    in.volumes_opt[vi].files.kix_path.c_str());
                kpx_cod[vi].close();
                kix_cod[vi].close();
                return {};
            }
            if (!kpx_opt[vi].open(in.volumes_opt[vi].files.kpx_path)) {
                if (logger) logger->error("Cannot open %s",
                    in.volumes_opt[vi].files.kpx_path.c_str());
                kix_opt[vi].close();
                kpx_cod[vi].close();
                kix_cod[vi].close();
                return {};
            }
        }

        bundles[vi].kix = &kix_cod[vi];
        bundles[vi].kpx = &kpx_cod[vi];
        if (in.both_mode) {
            bundles[vi].kix_opt = &kix_opt[vi];
            bundles[vi].kpx_opt = &kpx_opt[vi];
        }

        apply_stage2a_madvise(kix_cod[vi], kpx_cod[vi]);
        if (in.both_mode) {
            apply_stage2a_madvise(kix_opt[vi], kpx_opt[vi]);
        }
        s2_cur_vols.push_back(static_cast<uint16_t>(vi));
        s2_cur_ext_jobs += ext_jobs_per_volume;
        if (s2_cur_ext_jobs >= group_ext_job_target) {
            run_stage2_batch_and_close();
        }
    }
    run_stage2_batch_and_close();

    if (logger) {
        logger->info("Stage 2A complete: %zu hit(s) collected (aggregated)",
                     total_stage2a_hits);
        logger->info("Stage 2B complete: %zu chain(s)", results.size());
    }

    // Stage 2 dedup over the parent-relative key.  Fragment
    // splitting can yield duplicate chains when a query hits the same
    // parent region through several adjacent fragments; collapse them
    // here, before the (orch_hit -> output_hit) boundary, so neither
    // Stage 3 alignment nor TSV writing wastes work on duplicates.
    {
        const size_t before = results.size();
        dedup_stage2_orchestrator_hits(results, in.ksx_per_volume);
        if (logger && before != results.size()) {
            logger->info("Stage 2 dedup: %zu chain(s) -> %zu after dedup",
                         before, results.size());
        }
    }

    // Parent-level top-N chain selection.  After overlap dedup collapses
    // duplicate per-fragment chains, cap the chains kept per
    // (query, parent[, strand]) group by chainscore.
    {
        const size_t before = results.size();
        select_parent_topn(results, in.config.stage2.max_nhit_per_subject,
                           in.config.stage2.max_nhit_per_subject_mode,
                           in.ksx_per_volume);
        if (logger && before != results.size()) {
            logger->info("Parent top-N: %zu chain(s) -> %zu after selection",
                         before, results.size());
        }
    }

    // Stage 2 in-total (L) cap: per query, keep the top-L chains by chainscore
    // (tie-inclusive), bounding the chains that enter Stage 3 alignment.
    {
        const size_t before = results.size();
        apply_in_total_stage2(results, in.config.stage2.max_nhit_in_total);
        if (logger && before != results.size()) {
            logger->info("Stage 2 in-total: %zu chain(s) -> %zu after cap",
                         before, results.size());
        }
    }
    return results;
}

template std::vector<OrchestratorHit>
run_search<uint16_t>(const RunSearchInputs<uint16_t>&);
template std::vector<OrchestratorHit>
run_search<uint32_t>(const RunSearchInputs<uint32_t>&);

}  // namespace ikafssn
