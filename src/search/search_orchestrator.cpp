#include "search/search_orchestrator.hpp"

#include "index/kix_reader.hpp"
#include "index/kpx_reader.hpp"
#include "index/ksx_reader.hpp"
#include "search/parallel_search.hpp"
#include "search/oid_filter.hpp"
#include "search/result_dedup.hpp"
#include "search/stage1_filter.hpp"
#include "util/common_init.hpp"
#include "util/logger.hpp"

#include <algorithm>
#include <cstdio>
#include <iterator>
#include <utility>

#include <tbb/enumerable_thread_specific.h>
#include <tbb/task_arena.h>

namespace ikafssn {

namespace {

enum class Stage1Tier4 : uint8_t {
    kTier4 = 4,  // multi-volume WILLNEED
    kTier1 = 1,  // single-volume RANDOM
};

enum class Stage2ATier : uint8_t {
    kTier4 = 4,  // multi-volume WILLNEED for both kix and kpx
    kTier2 = 2,  // single-volume kpx WILLNEED + kix RANDOM
    kTier3 = 3,  // single-volume kix WILLNEED + kpx RANDOM
    kTier1 = 1,  // single-volume both RANDOM
};

struct Stage1Batch {
    std::vector<uint16_t> vols;
    Stage1Tier4 tier = Stage1Tier4::kTier4;
};

struct Stage2Batch {
    std::vector<uint16_t> vols;
    Stage2ATier tier = Stage2ATier::kTier4;
};

uint64_t kix_cost(const VolumeMeta& v_cod, const VolumeMeta* v_opt) {
    uint64_t c = v_cod.kix_full_size;
    if (v_opt) c += v_opt->kix_full_size;
    return c;
}

uint64_t kpx_cost(const VolumeMeta& v_cod, const VolumeMeta* v_opt) {
    uint64_t c = v_cod.kpx_full_size;
    if (v_opt) c += v_opt->kpx_full_size;
    return c;
}

std::vector<Stage1Batch>
plan_stage1_batches(const std::vector<VolumeMeta>& vols_cod,
                    const std::vector<VolumeMeta>& vols_opt,
                    uint64_t posting_budget,
                    bool both_mode) {
    std::vector<Stage1Batch> batches;
    Stage1Batch cur;
    cur.tier = Stage1Tier4::kTier4;
    uint64_t cur_size = 0;

    auto flush = [&]() {
        if (!cur.vols.empty()) {
            batches.push_back(std::move(cur));
            cur = Stage1Batch{};
            cur.tier = Stage1Tier4::kTier4;
            cur_size = 0;
        }
    };

    for (size_t vi = 0; vi < vols_cod.size(); ++vi) {
        const VolumeMeta* opt = both_mode ? &vols_opt[vi] : nullptr;
        uint64_t c = kix_cost(vols_cod[vi], opt);
        if (c <= posting_budget) {
            if (cur_size + c > posting_budget && !cur.vols.empty()) {
                flush();
            }
            cur.vols.push_back(static_cast<uint16_t>(vi));
            cur_size += c;
        } else {
            flush();
            Stage1Batch single;
            single.tier = Stage1Tier4::kTier1;
            single.vols.push_back(static_cast<uint16_t>(vi));
            batches.push_back(std::move(single));
        }
    }
    flush();
    return batches;
}

std::vector<Stage2Batch>
plan_stage2_batches(const std::vector<VolumeMeta>& vols_cod,
                     const std::vector<VolumeMeta>& vols_opt,
                     uint64_t posting_budget,
                     bool both_mode) {
    std::vector<Stage2Batch> batches;
    Stage2Batch cur;
    cur.tier = Stage2ATier::kTier4;
    uint64_t cur_size = 0;

    auto flush = [&]() {
        if (!cur.vols.empty()) {
            batches.push_back(std::move(cur));
            cur = Stage2Batch{};
            cur.tier = Stage2ATier::kTier4;
            cur_size = 0;
        }
    };

    for (size_t vi = 0; vi < vols_cod.size(); ++vi) {
        const VolumeMeta* opt = both_mode ? &vols_opt[vi] : nullptr;
        uint64_t c_full = kix_cost(vols_cod[vi], opt) + kpx_cost(vols_cod[vi], opt);
        if (c_full <= posting_budget) {
            if (cur_size + c_full > posting_budget && !cur.vols.empty()) {
                flush();
            }
            cur.vols.push_back(static_cast<uint16_t>(vi));
            cur_size += c_full;
        } else {
            flush();
            Stage2Batch single;
            single.vols.push_back(static_cast<uint16_t>(vi));

            uint64_t c_kpx = kpx_cost(vols_cod[vi], opt);
            uint64_t c_kix = kix_cost(vols_cod[vi], opt);
            if (c_kpx <= posting_budget) {
                single.tier = Stage2ATier::kTier2;
            } else if (c_kix <= posting_budget) {
                single.tier = Stage2ATier::kTier3;
            } else {
                single.tier = Stage2ATier::kTier1;
            }
            batches.push_back(std::move(single));
        }
    }
    flush();
    return batches;
}

void apply_stage1_madvise(KixReader& kix, Stage1Tier4 tier) {
    if (tier == Stage1Tier4::kTier4) {
        kix.apply_madvise_full(true);
    } else {
        kix.apply_madvise_posting_random();
    }
}

void release_stage1_madvise(KixReader& kix, Stage1Tier4 tier) {
    if (tier == Stage1Tier4::kTier4) {
        kix.apply_madvise_full(false);
    }
    // Tier 1: nothing to release; close() drops the mapping anyway.
}

void apply_stage2a_madvise(KixReader& kix, KpxReader& kpx, Stage2ATier tier) {
    switch (tier) {
        case Stage2ATier::kTier4:
            kix.apply_madvise_full(true);
            kpx.apply_madvise_full(true);
            break;
        case Stage2ATier::kTier2:
            kix.apply_madvise_posting_random();
            kpx.apply_madvise_full(true);
            break;
        case Stage2ATier::kTier3:
            kix.apply_madvise_full(true);
            kpx.apply_madvise_posting_random();
            break;
        case Stage2ATier::kTier1:
            kix.apply_madvise_posting_random();
            kpx.apply_madvise_posting_random();
            break;
    }
}

void release_stage2a_madvise(KixReader& kix, KpxReader& kpx, Stage2ATier tier) {
    switch (tier) {
        case Stage2ATier::kTier4:
            kix.apply_madvise_full(false);
            kpx.apply_madvise_full(false);
            break;
        case Stage2ATier::kTier2:
            kpx.apply_madvise_full(false);
            break;
        case Stage2ATier::kTier3:
            kix.apply_madvise_full(false);
            break;
        case Stage2ATier::kTier1:
            break;
    }
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
        buf.tier = in.tier;
        buf.ensure_capacity(in.max_num_seqs);
        return buf;
    };
    tbb::enumerable_thread_specific<Stage1Buffer> tls_bufs(make_tls_buf);

    // ----------------------------------------------------------------
    // Stage 1
    // ----------------------------------------------------------------
    auto s1_batches = plan_stage1_batches(in.volumes_cod, in.volumes_opt,
                                          in.posting_budget, in.both_mode);
    if (logger) {
        size_t n4 = 0, n1 = 0;
        for (auto& b : s1_batches) {
            if (b.tier == Stage1Tier4::kTier4) ++n4;
            else                                ++n1;
        }
        logger->info("Stage 1 batch plan: %zu batch(es) over %zu volume(s) "
                     "(tier4=%zu, tier1=%zu)",
                     s1_batches.size(), num_volumes, n4, n1);
    }

    for (auto& batch : s1_batches) {
        // Open kix readers for the batch.
        for (uint16_t vi : batch.vols) {
            if (!kix_cod[vi].open(in.volumes_cod[vi].files.kix_path)) {
                if (logger) logger->error("Cannot open %s",
                    in.volumes_cod[vi].files.kix_path.c_str());
                return {};
            }
            apply_stage1_madvise(kix_cod[vi], batch.tier);
            bundles[vi].kix = &kix_cod[vi];
            bundles[vi].kpx = nullptr;
            if (in.both_mode) {
                if (!kix_opt[vi].open(in.volumes_opt[vi].files.kix_path)) {
                    if (logger) logger->error("Cannot open %s",
                        in.volumes_opt[vi].files.kix_path.c_str());
                    return {};
                }
                apply_stage1_madvise(kix_opt[vi], batch.tier);
                bundles[vi].kix_opt = &kix_opt[vi];
                bundles[vi].kpx_opt = nullptr;
            }
        }

        auto idxs = indices_for_batch(ext_jobs, batch.vols);
        run_stage1_jobs<KmerInt>(idxs, ext_jobs, queries, bundles,
                                  in.k, in.config, in.both_mode,
                                  states, arena, tls_bufs);

        for (uint16_t vi : batch.vols) {
            release_stage1_madvise(kix_cod[vi], batch.tier);
            kix_cod[vi].close();
            bundles[vi].kix = nullptr;
            if (in.both_mode) {
                release_stage1_madvise(kix_opt[vi], batch.tier);
                kix_opt[vi].close();
                bundles[vi].kix_opt = nullptr;
            }
        }
    }

    if (logger) {
        size_t total_candidates = 0;
        for (auto& st : states) total_candidates += st.candidates.size();
        logger->info("Stage 1 complete: %zu candidate(s) across %zu ext_job(s)",
                     total_candidates, ext_jobs.size());
    }

    // ----------------------------------------------------------------
    // Mode 1 fast path
    // ----------------------------------------------------------------
    if (in.config.mode == 1) {
        std::vector<OrchestratorHit> results;
        for (size_t ji = 0; ji < ext_jobs.size(); ++ji) {
            const ExtJob& ej = ext_jobs[ji];
            JobState& st = states[ji];
            for (auto& cr : st.mode1_results) {
                OrchestratorHit oh;
                oh.query_idx = ej.qi;
                oh.volume_idx = ej.vi;
                oh.volume_index = volume_indices[ej.vi];
                oh.cr = std::move(cr);
                results.push_back(std::move(oh));
            }
            st.mode1_results.clear();
        }
        return results;
    }

    // ----------------------------------------------------------------
    // Stage 2 (Stage 2A + Stage 2B fused into one batch loop)
    // ----------------------------------------------------------------
    auto s2_batches = plan_stage2_batches(in.volumes_cod, in.volumes_opt,
                                          in.posting_budget, in.both_mode);
    if (logger) {
        size_t n4 = 0, n3 = 0, n2 = 0, n1 = 0;
        for (auto& b : s2_batches) {
            switch (b.tier) {
                case Stage2ATier::kTier4: ++n4; break;
                case Stage2ATier::kTier2: ++n2; break;
                case Stage2ATier::kTier3: ++n3; break;
                case Stage2ATier::kTier1: ++n1; break;
            }
        }
        logger->info("Stage 2 batch plan: %zu batch(es) over %zu volume(s) "
                     "(tier4=%zu, tier2=%zu, tier3=%zu, tier1=%zu)",
                     s2_batches.size(), num_volumes, n4, n2, n3, n1);
    }

    // Stage 2A and Stage 2B share this batch loop: Stage 2B runs immediately
    // after Stage 2A for the same batch and the per-ext_job transient state
    // is freed before the next batch starts, so peak memory tracks
    // `posting_budget` rather than the total volume × query fan-out.
    std::vector<OrchestratorHit> results;
    size_t total_stage2a_hits = 0;

    for (auto& batch : s2_batches) {
        for (uint16_t vi : batch.vols) {
            if (!kix_cod[vi].open(in.volumes_cod[vi].files.kix_path)) {
                if (logger) logger->error("Cannot open %s",
                    in.volumes_cod[vi].files.kix_path.c_str());
                return {};
            }
            if (!kpx_cod[vi].open(in.volumes_cod[vi].files.kpx_path)) {
                if (logger) logger->error("Cannot open %s",
                    in.volumes_cod[vi].files.kpx_path.c_str());
                return {};
            }
            apply_stage2a_madvise(kix_cod[vi], kpx_cod[vi], batch.tier);
            bundles[vi].kix = &kix_cod[vi];
            bundles[vi].kpx = &kpx_cod[vi];
            if (in.both_mode) {
                if (!kix_opt[vi].open(in.volumes_opt[vi].files.kix_path)) {
                    if (logger) logger->error("Cannot open %s",
                        in.volumes_opt[vi].files.kix_path.c_str());
                    return {};
                }
                if (!kpx_opt[vi].open(in.volumes_opt[vi].files.kpx_path)) {
                    if (logger) logger->error("Cannot open %s",
                        in.volumes_opt[vi].files.kpx_path.c_str());
                    return {};
                }
                apply_stage2a_madvise(kix_opt[vi], kpx_opt[vi], batch.tier);
                bundles[vi].kix_opt = &kix_opt[vi];
                bundles[vi].kpx_opt = &kpx_opt[vi];
            }
        }

        auto idxs = indices_for_batch(ext_jobs, batch.vols);
        run_stage2a_jobs<KmerInt>(idxs, ext_jobs, queries, bundles,
                                   in.both_mode, states, arena);

        // Close .kix/.kpx readers before Stage 2B — Stage 2B reads only
        // JobState and does not need the posting files.
        for (uint16_t vi : batch.vols) {
            release_stage2a_madvise(kix_cod[vi], kpx_cod[vi], batch.tier);
            kpx_cod[vi].close();
            kix_cod[vi].close();
            bundles[vi].kix = nullptr;
            bundles[vi].kpx = nullptr;
            if (in.both_mode) {
                release_stage2a_madvise(kix_opt[vi], kpx_opt[vi], batch.tier);
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
            logger->info("Stage 2A batch (%zu vol(s), %zu ext_job(s)): %zu hit(s)",
                         batch.vols.size(), idxs.size(), batch_hits);
        }

        size_t chains_before = results.size();
        run_stage2b_jobs(idxs, ext_jobs, states, volume_indices, arena, results);
        if (logger) {
            logger->info("Stage 2B batch: %zu chain(s)",
                         results.size() - chains_before);
        }

        // Release this batch's per-ext_job transient state to bound peak
        // memory.  Move-assign from {} (instead of clear()) so the bucket
        // arrays / vector capacity are actually returned to the allocator.
        for (size_t ji : idxs) {
            JobState& st = states[ji];
            st.hits_per_seq  = {};
            st.stage1_scores = {};
            std::vector<SeqId>().swap(st.sorted_candidate_sids);
            std::vector<Stage1Candidate>().swap(st.candidates);
        }
    }

    if (logger) {
        logger->info("Stage 2A complete: %zu hit(s) collected (aggregated)",
                     total_stage2a_hits);
        logger->info("Stage 2B complete: %zu chain(s)", results.size());
    }

    // v10 (Phase 3): Stage 2 dedup over the parent-OID-relative key.
    // Fragment splitting can yield duplicate chains when a query hits the
    // same parent region through several adjacent fragments; collapse them
    // here, before the (orch_hit -> output_hit) boundary, so neither Stage 3
    // alignment nor TSV writing wastes work on duplicates.
    {
        const size_t before = results.size();
        dedup_stage2_orchestrator_hits(results, in.ksx_per_volume);
        if (logger && before != results.size()) {
            logger->info("Stage 2 dedup: %zu chain(s) -> %zu after dedup",
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
