// bench_stage1 — measure Stage 1 batch SIMD update kernel throughput per tier.
//
// The Stage 1 hot loop calls flush_batch_simd<Tier> with sid batches; this
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
        flush_batch_simd<Stage1Tier::T32>(
            sids.data(), batch_size, /*q_pos=*/42u,
            static_cast<void*>(scores.data()),
            static_cast<void*>(last_pos.data()),
            dirty);

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
        static_cast<int>(SimdCap::Scalar),
        static_cast<int>(SimdCap::SSE42),
        static_cast<int>(SimdCap::AVX2),
        static_cast<int>(SimdCap::AVX512BW),
        static_cast<int>(SimdCap::AVX512VBMI),
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

int main(int argc, char** argv) {
    init_simd_dispatch(nullptr);
    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
    return 0;
}
