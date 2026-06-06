#include "search/search_orchestrator.hpp"

#include "index/kix_reader.hpp"
#include "index/kpx_reader.hpp"
#include "index/ksx_reader.hpp"
#include "search/decode_cache.hpp"
#include "search/parallel_search.hpp"
#include "search/oid_filter.hpp"
#include "search/query_preprocessor.hpp"
#include "search/result_dedup.hpp"
#include "search/stage1_filter.hpp"
#include "util/common_init.hpp"
#include "util/logger.hpp"

#include <algorithm>
#include <cstdio>
#include <iterator>
#include <unordered_set>
#include <utility>

#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_scan.h>
#include <tbb/task_arena.h>

namespace ikafssn {

namespace {

// Stage 1 posting-body madvise mode.
enum class Stage1Madvise : uint8_t {
    kBatchWillneed,  // multi-volume batch: WILLNEED the needed ranges
    kSoloRandom,     // single over-budget volume: posting body RANDOM
};

// Stage 2A madvise mode: which of kix / kpx get WILLNEED vs RANDOM.
enum class Stage2AMadvise : uint8_t {
    kBothWillneed,
    kKpxWillneedKixRandom,
    kKixWillneedKpxRandom,
    kBothRandom,
};

inline const char* madvise_name(Stage1Madvise m) {
    return m == Stage1Madvise::kBatchWillneed ? "batch-willneed" : "solo-random";
}
inline const char* madvise_name(Stage2AMadvise m) {
    switch (m) {
        case Stage2AMadvise::kBothWillneed:         return "both-willneed";
        case Stage2AMadvise::kKpxWillneedKixRandom: return "kpx-willneed/kix-random";
        case Stage2AMadvise::kKixWillneedKpxRandom: return "kix-willneed/kpx-random";
        case Stage2AMadvise::kBothRandom:           return "both-random";
    }
    return "?";
}

// (offset, length) into a posting file.
using ByteRange = std::pair<uint64_t, uint64_t>;

// Page-size gap threshold for ByteRange coalescing. Smaller gaps merge so
// that one MADV_WILLNEED call covers contiguous-or-near-contiguous posting
// lists, amortising the syscall over the prefetch window.
constexpr uint64_t kCoalesceGap = 4096;

void coalesce_ranges(std::vector<ByteRange>& ranges) {
    if (ranges.empty()) return;
    std::sort(ranges.begin(), ranges.end());
    size_t w = 0;
    for (size_t r = 0; r < ranges.size(); ++r) {
        if (w == 0) {
            ranges[w++] = ranges[r];
            continue;
        }
        ByteRange& last = ranges[w - 1];
        uint64_t last_end = last.first + last.second;
        if (ranges[r].first <= last_end + kCoalesceGap) {
            uint64_t new_end = std::max(last_end,
                                        ranges[r].first + ranges[r].second);
            last.second = new_end - last.first;
        } else {
            ranges[w++] = ranges[r];
        }
    }
    ranges.resize(w);
}

uint64_t total_range_bytes(const std::vector<ByteRange>& ranges) {
    uint64_t total = 0;
    for (auto& r : ranges) total += r.second;
    return total;
}

// Collect the set of unique k-mers each volume must probe in Stage 1 and
// Stage 2A.  `which` selects fwd / rc / both based on config.strand.
// `secondary` selects the optimal-side kmer stream (for both-mode) when
// true. Returns a sorted vector of unique k-mer values.
template <typename KmerInt>
std::vector<KmerInt> collect_unique_kmers(
    const std::vector<QueryBundle<KmerInt>>& queries,
    const std::vector<uint8_t>* skip_reason,
    int8_t strand,
    bool secondary) {
    // Mirror build_ext_jobs's strand handling: 2 = both, 1 = fwd only,
    // everything else (including the -1 / 0 server-default) = rc only.
    const bool want_fwd = (strand == 1 || strand == 2);
    const bool want_rc  = !want_fwd || (strand == 2);
    std::unordered_set<KmerInt> seen;
    auto take_kmers = [&](const QueryKmerData<KmerInt>* qd) {
        if (!qd) return;
        if (want_fwd) {
            for (auto k : qd->fwd_kmer_values) seen.insert(k);
        }
        if (want_rc) {
            for (auto k : qd->rc_kmer_values) seen.insert(k);
        }
    };
    for (size_t qi = 0; qi < queries.size(); ++qi) {
        if (skip_reason && qi < skip_reason->size() && (*skip_reason)[qi] != 0)
            continue;
        const auto* qd = secondary ? queries[qi].qdata_secondary
                                   : queries[qi].qdata_primary;
        take_kmers(qd);
    }
    std::vector<KmerInt> out(seen.begin(), seen.end());
    std::sort(out.begin(), out.end());
    return out;
}

// Walk an already-open kix volume's dictionary for every unique k-mer and
// build the coalesced (offset, length) posting-list ranges. Operates on
// the caller's open KixReader so the dict pages stay hot for the
// subsequent madvise + search on the same volume.
template <typename KmerInt>
std::vector<ByteRange> compute_kix_ranges(
    const KixReader& kix,
    const std::vector<KmerInt>& unique_kmers) {
    std::vector<ByteRange> ranges;
    ranges.reserve(unique_kmers.size());
    for (KmerInt k : unique_kmers) {
        uint64_t off, len;
        kix.posting_list_range(static_cast<uint32_t>(k), off, len);
        if (len == 0) continue;
        ranges.emplace_back(off, len);
    }
    coalesce_ranges(ranges);
    return ranges;
}

template <typename KmerInt>
std::vector<ByteRange> compute_kpx_ranges(
    const KpxReader& kpx,
    const std::vector<KmerInt>& unique_kmers) {
    std::vector<ByteRange> ranges;
    ranges.reserve(unique_kmers.size());
    for (KmerInt k : unique_kmers) {
        uint64_t lo, hi;
        kpx.pos_offset_range(static_cast<uint32_t>(k), lo, hi);
        if (hi <= lo) continue;
        ranges.emplace_back(lo, hi - lo);
    }
    coalesce_ranges(ranges);
    return ranges;
}

// Stage 1 madvise: dict is always pinned (WILLNEED + HUGEPAGE); the posting
// body is either prefetched as coalesced ranges (kBatchWillneed — the batch
// fits the budget) or marked RANDOM (kSoloRandom — one volume's needed bytes
// exceed the budget). The ranges are the exact (offset, length) pairs the
// batch will probe, so the kernel's readahead window tracks real demand
// instead of the full posting file.
void apply_stage1_madvise(KixReader& kix, Stage1Madvise mode,
                          const std::vector<ByteRange>& ranges) {
    if (mode == Stage1Madvise::kBatchWillneed) {
        kix.apply_madvise_dict_only();
        kix.apply_madvise_posting_ranges(ranges);
    } else {
        kix.apply_madvise_posting_random();
    }
}

void apply_stage2a_madvise(KixReader& kix, KpxReader& kpx, Stage2AMadvise mode,
                           const std::vector<ByteRange>& kix_ranges,
                           const std::vector<ByteRange>& kpx_ranges) {
    switch (mode) {
        case Stage2AMadvise::kBothWillneed:
            kix.apply_madvise_dict_only();
            kpx.apply_madvise_dict_only();
            kix.apply_madvise_posting_ranges(kix_ranges);
            kpx.apply_madvise_posting_ranges(kpx_ranges);
            break;
        case Stage2AMadvise::kKpxWillneedKixRandom:
            kix.apply_madvise_posting_random();
            kpx.apply_madvise_dict_only();
            kpx.apply_madvise_posting_ranges(kpx_ranges);
            break;
        case Stage2AMadvise::kKixWillneedKpxRandom:
            kix.apply_madvise_dict_only();
            kix.apply_madvise_posting_ranges(kix_ranges);
            kpx.apply_madvise_posting_random();
            break;
        case Stage2AMadvise::kBothRandom:
            kix.apply_madvise_posting_random();
            kpx.apply_madvise_posting_random();
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

    // Per-volume Stage 1 decode caches (one slot per volume; addresses stay
    // stable so set_decode_cache pointers remain valid).  A volume's cache is
    // materialised in the open loop when it joins a batch and released right
    // after the batch runs, so resident decode heap is bounded by the batch's
    // charged decode bytes.  s1_cacheable[vi] records whether vi's cache fit
    // the budget and was materialised.
    std::vector<DecodedKmerCache> caches_cod(num_volumes), caches_opt;
    if (in.both_mode) caches_opt.resize(num_volumes);
    std::vector<uint8_t> s1_cacheable(num_volumes, 0);

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

    // ----------------------------------------------------------------
    // Pre-compute the set of unique k-mers each batch must probe.  The
    // Stage 1 / Stage 2 batch loops walk every volume's dictionary
    // incrementally and decide on the fly whether the next volume fits
    // in the current batch's posting_budget.
    // ----------------------------------------------------------------
    auto unique_cod = collect_unique_kmers<KmerInt>(
        queries, in.query_skip_reason, in.config.strand, /*secondary=*/false);
    std::vector<KmerInt> unique_opt;
    if (in.both_mode) {
        unique_opt = collect_unique_kmers<KmerInt>(
            queries, in.query_skip_reason, in.config.strand, /*secondary=*/true);
    }
    if (logger) {
        logger->info("Search planning: %zu unique kmers (cod=%zu, opt=%zu) over %zu volume(s); posting_budget=%llu",
                     unique_cod.size() + unique_opt.size(),
                     unique_cod.size(), unique_opt.size(), num_volumes,
                     static_cast<unsigned long long>(in.posting_budget));
    }

    // Mode 1 fast path uses this batch-local accumulator: each Stage 1
    // batch flushes its mode1_results into it and releases the per-ji
    // capacity, so the orchestrator's peak heap stays bounded by one
    // batch's worth of candidates instead of all batches combined.
    std::vector<OrchestratorHit> mode1_results;

    // ----------------------------------------------------------------
    // Stage 1: incremental open + WILLNEED + accumulate.  Each
    // volume's dict is walked once for the cost and prefetch decision;
    // the same KixReader stays open for the search that follows.
    // ----------------------------------------------------------------
    std::vector<uint16_t> s1_cur_vols;
    uint64_t              s1_cur_bytes = 0;
    size_t                s1_batch_seq = 0;

    auto run_stage1_batch_and_close = [&](Stage1Madvise mode) {
        if (s1_cur_vols.empty()) return;
        auto idxs = indices_for_batch(ext_jobs, s1_cur_vols);

        // Bind the decode cache of every cacheable volume in the batch (the
        // caches were materialised in the open loop).  set_decode_cache is a
        // per-KixReader handle, so a batch may hold a mix of cached and
        // uncached volumes; the Stage 1 hot loop reads kix.decode_cache() per
        // volume and falls back to on-the-fly decode where unbound.  All batch
        // volumes then run in one parallel pass, so cross-volume parallelism
        // scales with how many volumes the memory budget admits.
        for (uint16_t vi : s1_cur_vols) {
            if (!s1_cacheable[vi]) continue;
            kix_cod[vi].set_decode_cache(&caches_cod[vi]);
            if (in.both_mode) kix_opt[vi].set_decode_cache(&caches_opt[vi]);
        }

        run_stage1_jobs<KmerInt>(idxs, ext_jobs, queries, bundles,
                                  in.k, in.config, in.both_mode,
                                  states, arena, tls_bufs);

        // Unbind and free the batch's decode caches now that the pass is done,
        // bounding resident decode heap to one batch's charged bytes.
        for (uint16_t vi : s1_cur_vols) {
            if (!s1_cacheable[vi]) continue;
            kix_cod[vi].set_decode_cache(nullptr);
            caches_cod[vi].release();
            if (in.both_mode) {
                kix_opt[vi].set_decode_cache(nullptr);
                caches_opt[vi].release();
            }
        }

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
            logger->info("Stage 1 batch %zu: %zu vol(s), %zu ext_job(s), madvise=%s, charged_bytes=%llu (WILLNEED + decode heap)",
                         s1_batch_seq, s1_cur_vols.size(), idxs.size(),
                         madvise_name(mode),
                         static_cast<unsigned long long>(s1_cur_bytes));
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
        s1_cur_bytes = 0;
    };

    for (size_t vi = 0; vi < num_volumes; ++vi) {
        // Open the volume's kix (cod + opt when in both-mode) so the
        // cost computation and the subsequent madvise / search all share
        // the same open mapping — no double open per volume.
        if (!kix_cod[vi].open(in.volumes_cod[vi].files.kix_path)) {
            if (logger) logger->error("Cannot open %s",
                in.volumes_cod[vi].files.kix_path.c_str());
            return {};
        }
        std::vector<ByteRange> ranges_cod =
            compute_kix_ranges<KmerInt>(kix_cod[vi], unique_cod);
        uint64_t vi_bytes = kix_cod[vi].dict_size() + total_range_bytes(ranges_cod);

        std::vector<ByteRange> ranges_opt;
        if (in.both_mode) {
            if (!kix_opt[vi].open(in.volumes_opt[vi].files.kix_path)) {
                if (logger) logger->error("Cannot open %s",
                    in.volumes_opt[vi].files.kix_path.c_str());
                kix_cod[vi].close();
                return {};
            }
            ranges_opt = compute_kix_ranges<KmerInt>(kix_opt[vi], unique_opt);
            vi_bytes += kix_opt[vi].dict_size() + total_range_bytes(ranges_opt);
        }

        bundles[vi].kix = &kix_cod[vi];
        bundles[vi].kpx = nullptr;
        if (in.both_mode) {
            bundles[vi].kix_opt = &kix_opt[vi];
            bundles[vi].kpx_opt = nullptr;
        }

        if (vi_bytes > in.posting_budget) {
            // Over-budget volume: the WILLNEED bytes alone exceed the budget,
            // so drain the current batch, then run this volume on its own with
            // MADV_RANDOM on the posting body and no decode cache (it is
            // already memory-pressured).
            run_stage1_batch_and_close(Stage1Madvise::kBatchWillneed);

            apply_stage1_madvise(kix_cod[vi], Stage1Madvise::kSoloRandom, ranges_cod);
            if (in.both_mode) {
                apply_stage1_madvise(kix_opt[vi], Stage1Madvise::kSoloRandom, ranges_opt);
            }
            s1_cacheable[vi] = 0;
            s1_cur_vols.push_back(static_cast<uint16_t>(vi));
            s1_cur_bytes = vi_bytes;
            run_stage1_batch_and_close(Stage1Madvise::kSoloRandom);
        } else {
            // Apply the volume's WILLNEED first so the decode cache fill below
            // (header sizing + parallel decode) rides the prefetch.
            apply_stage1_madvise(kix_cod[vi], Stage1Madvise::kBatchWillneed, ranges_cod);
            if (in.both_mode) {
                apply_stage1_madvise(kix_opt[vi], Stage1Madvise::kBatchWillneed, ranges_opt);
            }

            // Build the decode cache and charge its heap against the same
            // (soft) budget as the WILLNEED page cache.  fill() returns false
            // for an over-budget or corrupt volume, which then runs uncached
            // (on-the-fly decode) but is still WILLNEED-batched.  In both-mode
            // the opt cache must fit alongside the cod cache, so it is filled
            // with the residual budget to keep the pair within one budget.
            uint64_t decoded_vi = 0;
            bool cacheable = caches_cod[vi].fill(kix_cod[vi], unique_cod,
                                                 arena, in.posting_budget);
            if (cacheable) decoded_vi = caches_cod[vi].decoded_bytes();
            if (cacheable && in.both_mode) {
                size_t remaining = (in.posting_budget > decoded_vi)
                    ? in.posting_budget - decoded_vi : 0;
                if (caches_opt[vi].fill(kix_opt[vi], unique_opt, arena, remaining)) {
                    decoded_vi += caches_opt[vi].decoded_bytes();
                } else {
                    caches_cod[vi].release();
                    cacheable = false;
                    decoded_vi = 0;
                }
            }
            s1_cacheable[vi] = cacheable ? 1 : 0;
            if (logger) {
                logger->info("Stage 1 vol %zu: decode cache %s (%llu decoded byte(s))",
                             vi, cacheable ? "cached" : "uncached",
                             static_cast<unsigned long long>(decoded_vi));
            }

            s1_cur_vols.push_back(static_cast<uint16_t>(vi));
            // Charge WILLNEED page cache + decode heap against the soft budget.
            // Add-then-check: the volume that crosses joins this batch (the
            // tolerated one-volume overshoot of the soft -memory_limit); more
            // volumes per batch as the budget grows raises cross-volume
            // parallelism.  Caches are released when the batch closes.
            s1_cur_bytes += vi_bytes + decoded_vi;
            if (s1_cur_bytes >= in.posting_budget) {
                run_stage1_batch_and_close(Stage1Madvise::kBatchWillneed);
            }
        }
    }
    run_stage1_batch_and_close(Stage1Madvise::kBatchWillneed);

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
        return mode1_results;
    }

    // ----------------------------------------------------------------
    // Stage 2 (Stage 2A + Stage 2B fused into one incremental loop)
    //
    // Stage 2A and Stage 2B share this loop: Stage 2B runs immediately
    // after Stage 2A for the same batch and the per-ext_job transient
    // state is freed before the next batch starts, so peak memory tracks
    // posting_budget rather than the total volume × query fan-out.
    // ----------------------------------------------------------------
    std::vector<OrchestratorHit> results;
    size_t total_stage2a_hits = 0;
    std::vector<uint16_t> s2_cur_vols;
    uint64_t              s2_cur_bytes = 0;
    size_t                s2_batch_seq = 0;

    auto run_stage2_batch_and_close = [&](Stage2AMadvise mode) {
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
            logger->info("Stage 2A batch %zu: %zu vol(s), %zu ext_job(s), madvise=%s, needed_bytes=%llu, %zu hit(s)",
                         s2_batch_seq, s2_cur_vols.size(), idxs.size(),
                         madvise_name(mode),
                         static_cast<unsigned long long>(s2_cur_bytes),
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
        s2_cur_bytes = 0;
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
        std::vector<ByteRange> kix_ranges_cod =
            compute_kix_ranges<KmerInt>(kix_cod[vi], unique_cod);
        std::vector<ByteRange> kpx_ranges_cod =
            compute_kpx_ranges<KmerInt>(kpx_cod[vi], unique_cod);
        uint64_t kix_bytes = kix_cod[vi].dict_size() + total_range_bytes(kix_ranges_cod);
        uint64_t kpx_bytes = kpx_cod[vi].dict_size() + total_range_bytes(kpx_ranges_cod);

        std::vector<ByteRange> kix_ranges_opt;
        std::vector<ByteRange> kpx_ranges_opt;
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
            kix_ranges_opt = compute_kix_ranges<KmerInt>(kix_opt[vi], unique_opt);
            kpx_ranges_opt = compute_kpx_ranges<KmerInt>(kpx_opt[vi], unique_opt);
            kix_bytes += kix_opt[vi].dict_size() + total_range_bytes(kix_ranges_opt);
            kpx_bytes += kpx_opt[vi].dict_size() + total_range_bytes(kpx_ranges_opt);
        }

        bundles[vi].kix = &kix_cod[vi];
        bundles[vi].kpx = &kpx_cod[vi];
        if (in.both_mode) {
            bundles[vi].kix_opt = &kix_opt[vi];
            bundles[vi].kpx_opt = &kpx_opt[vi];
        }

        const uint64_t vi_bytes = kix_bytes + kpx_bytes;

        if (vi_bytes > in.posting_budget) {
            // Over-budget volume: drain the accumulated batch first, then pick
            // the tightest madvise mode that fits the budget.
            run_stage2_batch_and_close(Stage2AMadvise::kBothWillneed);

            Stage2AMadvise solo_mode;
            if (kpx_bytes <= in.posting_budget) {
                solo_mode = Stage2AMadvise::kKpxWillneedKixRandom;
            } else if (kix_bytes <= in.posting_budget) {
                solo_mode = Stage2AMadvise::kKixWillneedKpxRandom;
            } else {
                solo_mode = Stage2AMadvise::kBothRandom;
            }

            apply_stage2a_madvise(kix_cod[vi], kpx_cod[vi], solo_mode,
                                  kix_ranges_cod, kpx_ranges_cod);
            if (in.both_mode) {
                apply_stage2a_madvise(kix_opt[vi], kpx_opt[vi], solo_mode,
                                      kix_ranges_opt, kpx_ranges_opt);
            }
            s2_cur_vols.push_back(static_cast<uint16_t>(vi));
            s2_cur_bytes = vi_bytes;
            run_stage2_batch_and_close(solo_mode);
        } else {
            apply_stage2a_madvise(kix_cod[vi], kpx_cod[vi], Stage2AMadvise::kBothWillneed,
                                  kix_ranges_cod, kpx_ranges_cod);
            if (in.both_mode) {
                apply_stage2a_madvise(kix_opt[vi], kpx_opt[vi], Stage2AMadvise::kBothWillneed,
                                      kix_ranges_opt, kpx_ranges_opt);
            }
            s2_cur_vols.push_back(static_cast<uint16_t>(vi));
            s2_cur_bytes += vi_bytes;
            if (s2_cur_bytes >= in.posting_budget) {
                run_stage2_batch_and_close(Stage2AMadvise::kBothWillneed);
            }
        }
    }
    run_stage2_batch_and_close(Stage2AMadvise::kBothWillneed);

    if (logger) {
        logger->info("Stage 2A complete: %zu hit(s) collected (aggregated)",
                     total_stage2a_hits);
        logger->info("Stage 2B complete: %zu chain(s)", results.size());
    }

    // Stage 2 dedup over the parent-OID-relative key.  Fragment
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
    return results;
}

template std::vector<OrchestratorHit>
run_search<uint16_t>(const RunSearchInputs<uint16_t>&);
template std::vector<OrchestratorHit>
run_search<uint32_t>(const RunSearchInputs<uint32_t>&);

}  // namespace ikafssn
