// bench_parallel_sort — measure index_builder's posting-buffer sort per tier.
//
// Compares scalar (tbb::parallel_sort) against the SIMD path (x86-simd-sort
// keyvalue_qsort + local pos sort). Each tier is forced via force_simd_cap()
// inside the bench fixture.

#include "bench/bench_common.hpp"
#include "index/parallel_sort_dispatch.hpp"
#include "util/simd_dispatch.hpp"

#include <benchmark/benchmark.h>

#include <cstdint>
#include <random>
#include <vector>

using namespace ikafssn;

namespace {

void BM_ParallelSortTempEntries(benchmark::State& state) {
    SimdCap tier = static_cast<SimdCap>(state.range(0));
    const std::size_t n = static_cast<std::size_t>(state.range(1));

    bench::apply_force_tier(state, tier);
    bench::label_tier(state, tier);

    std::mt19937 rng(0xCAFEBABE);
    std::vector<TempEntry> orig;
    orig.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        orig.push_back({static_cast<uint32_t>(rng() & 0xFFFFFFu),
                        static_cast<uint32_t>(rng() & 0xFFFFu),
                        static_cast<uint32_t>(rng() & 0xFFFFu)});
    }

    std::vector<TempEntry> buf;
    buf.reserve(n);

    for (auto _ : state) {
        state.PauseTiming();
        buf = orig;
        state.ResumeTiming();
        parallel_sort_temp_entries(buf);
        benchmark::DoNotOptimize(buf.data());
    }
    state.SetItemsProcessed(static_cast<int64_t>(state.iterations()) * n);
    state.SetBytesProcessed(
        static_cast<int64_t>(state.iterations()) * n * sizeof(TempEntry));
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
    static constexpr int64_t kSizes[] = {
        1 << 10, 1 << 14, 1 << 18, 1 << 22  // 1K, 16K, 256K, 4M
    };
    for (int t : kTiers) {
        for (int64_t n : kSizes) b->Args({t, n});
    }
}

} // namespace

BENCHMARK(BM_ParallelSortTempEntries)->Apply(RegisterArgs)->UseRealTime();

int main(int argc, char** argv) {
    init_simd_dispatch(nullptr);
    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
    return 0;
}
