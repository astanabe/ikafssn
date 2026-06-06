// bench_stage1 — measure Stage 1 batch SIMD update kernel throughput per tier.
//
// The Stage 1 hot loop calls flush_batch_simd<Width> with sid batches; this
// benchmark calls the kernel directly with synthetic batches so per-tier
// dispatch overhead, gather/scatter cost, and the conflict-detection scalar
// fallback can be measured in isolation.

#include "bench/bench_common.hpp"
#include "search/stage1_filter.hpp"
#include "search/stage1_filter_simd.hpp"
#include "util/simd_dispatch.hpp"

#include <benchmark/benchmark.h>

#include <cstdint>
#include <cstring>
#include <limits>
#include <random>
#include <vector>

using namespace ikafssn;

namespace {

void BM_Stage1FlushBatchT32(benchmark::State& state) {
    SimdCap tier = static_cast<SimdCap>(state.range(0));
    const int    batch_size = static_cast<int>(state.range(1));
    const size_t num_seqs   = static_cast<size_t>(state.range(2));

    bench::apply_force_tier(state, tier);
    bench::label_tier(state, tier);

    std::vector<uint32_t> scores(num_seqs, 0);
    std::vector<uint32_t> last_pos(num_seqs, std::numeric_limits<uint32_t>::max());
    std::vector<SeqId>    sids(static_cast<size_t>(batch_size));

    std::mt19937 rng(0xC0FFEE);
    for (auto& s : sids) s = static_cast<SeqId>(rng() % num_seqs);

    std::vector<uint32_t> dirty;
    dirty.reserve(static_cast<size_t>(batch_size));

    for (auto _ : state) {
        flush_batch_simd<Stage1Width::T32>(
            sids.data(), batch_size, /*q_pos=*/42u,
            static_cast<void*>(scores.data()),
            static_cast<void*>(last_pos.data()),
            dirty, /*remaining=*/0u, /*cutoff_T=*/0u);

        // Reset state without skewing the timed region.
        state.PauseTiming();
        std::memset(scores.data(), 0, scores.size() * sizeof(uint32_t));
        std::memset(last_pos.data(), 0xFF, last_pos.size() * sizeof(uint32_t));
        dirty.clear();
        state.ResumeTiming();
    }
    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) * batch_size);
}

void RegisterArgs(benchmark::internal::Benchmark* b) {
    static constexpr int kTiers[] = {
        static_cast<int>(SimdCap::SSE42),
        static_cast<int>(SimdCap::AVX2),
        static_cast<int>(SimdCap::AVX512BW),
        static_cast<int>(SimdCap::AVX512VBMI2),
    };
    static constexpr int kBatchSizes[]   = {16, 32, 64};
    static constexpr int64_t kNumSeqs[]  = {4096, 65536, 1 << 20};
    for (int t : kTiers) {
        for (int bs : kBatchSizes) {
            for (int64_t ns : kNumSeqs) b->Args({t, bs, ns});
        }
    }
}

} // namespace

BENCHMARK(BM_Stage1FlushBatchT32)->Apply(RegisterArgs)->UseRealTime();

// Measure clear_dirty_typed cost across dirty ratios. The kernel takes the
// per-index reset path below the threshold (dirty.size()*8 <= capacity) and
// the bulk reset_all path above it; both are exercised here.
void BM_Stage1ClearDirtyT32(benchmark::State& state) {
    const size_t num_seqs    = static_cast<size_t>(state.range(0));
    const size_t num_dirty   = static_cast<size_t>(state.range(1));

    Stage1Buffer buf;
    buf.width = Stage1Width::T32;
    buf.ensure_capacity(static_cast<uint32_t>(num_seqs));

    std::mt19937 rng(0xBADF00D);
    std::vector<uint32_t> sample;
    sample.reserve(num_dirty);
    for (size_t i = 0; i < num_dirty; i++) {
        sample.push_back(static_cast<uint32_t>(rng() % num_seqs));
    }

    auto* scores   = score_ptr<Stage1Width::T32>(buf);
    auto* last_pos = last_pos_ptr<Stage1Width::T32>(buf);

    for (auto _ : state) {
        state.PauseTiming();
        buf.dirty.clear();
        for (uint32_t sid : sample) {
            scores[sid]   = 1;
            last_pos[sid] = 7;
            buf.dirty.push_back(sid);
        }
        state.ResumeTiming();
        buf.clear_dirty_typed<Stage1Width::T32>();
    }
    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) * num_dirty);
}

void RegisterClearDirtyArgs(benchmark::internal::Benchmark* b) {
    static constexpr int64_t kNumSeqs[]  = {1 << 16, 1 << 20};
    // dirty fractions: 1/64, 1/16, 1/8, 1/2 — straddling the bulk-reset threshold.
    static constexpr int kDirtyDenoms[]  = {64, 16, 8, 2};
    for (int64_t ns : kNumSeqs) {
        for (int denom : kDirtyDenoms) {
            b->Args({ns, ns / denom});
        }
    }
}

BENCHMARK(BM_Stage1ClearDirtyT32)->Apply(RegisterClearDirtyArgs)->UseRealTime();

int main(int argc, char** argv) {
    init_simd_dispatch(nullptr);
    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
    return 0;
}
