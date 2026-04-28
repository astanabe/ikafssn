// Phase 3b — k-mer reverse-complement batch SIMD micro-benchmark.
//
// Measures kmer_revcomp_batch() at each available SIMD tier vs the scalar
// kmer_revcomp() loop, for representative k values and query lengths.

#include "core/kmer_encoding.hpp"
#include "core/kmer_revcomp_simd.hpp"
#include "util/simd_dispatch.hpp"
#include "bench_common.hpp"

#include <benchmark/benchmark.h>

#include <cstdint>
#include <random>
#include <type_traits>
#include <vector>

using namespace ikafssn;

// Scalar baseline: per-element kmer_revcomp() in a loop.
template <typename KmerInt>
static void BM_KmerRevcompScalarLoop(benchmark::State& state) {
    int k = static_cast<int>(state.range(0));
    std::size_t n = static_cast<std::size_t>(state.range(1));
    std::mt19937 rng(42);
    KmerInt mask = kmer_mask<KmerInt>(k);
    std::vector<KmerInt> in(n), out(n);
    for (auto& v : in) v = static_cast<KmerInt>(rng()) & mask;
    for (auto _ : state) {
        for (std::size_t i = 0; i < n; i++) out[i] = kmer_revcomp<KmerInt>(in[i], k);
        benchmark::DoNotOptimize(out.data());
    }
    state.SetItemsProcessed(int64_t(state.iterations()) * n);
}
BENCHMARK(BM_KmerRevcompScalarLoop<uint16_t>)->ArgsProduct({{6, 8}, {32, 256, 4096}});
BENCHMARK(BM_KmerRevcompScalarLoop<uint32_t>)->ArgsProduct({{11, 12}, {32, 256, 4096}});

template <typename KmerInt>
static void BM_KmerRevcompBatch(benchmark::State& state) {
    SimdCap tier = static_cast<SimdCap>(state.range(0));
    int k = static_cast<int>(state.range(1));
    std::size_t n = static_cast<std::size_t>(state.range(2));
    bench::apply_force_tier(state, tier);
    bench::label_tier(state, tier);

    std::mt19937 rng(42);
    KmerInt mask = kmer_mask<KmerInt>(k);
    std::vector<KmerInt> in(n), out(n);
    for (auto& v : in) v = static_cast<KmerInt>(rng()) & mask;
    for (auto _ : state) {
        kmer_revcomp_batch<KmerInt>(in.data(), out.data(), n, k);
        benchmark::DoNotOptimize(out.data());
    }
    state.SetItemsProcessed(int64_t(state.iterations()) * n);
}

// uint16_t (k <= 8)
BENCHMARK(BM_KmerRevcompBatch<uint16_t>)
    ->ArgsProduct({{
        (int64_t)SimdCap::SSE42,
        (int64_t)SimdCap::AVX2,
        (int64_t)SimdCap::AVX512BW,
        (int64_t)SimdCap::AVX512VBMI2,
    }, {6, 8}, {32, 256, 4096}});

// uint32_t (k = 9..16)
BENCHMARK(BM_KmerRevcompBatch<uint32_t>)
    ->ArgsProduct({{
        (int64_t)SimdCap::SSE42,
        (int64_t)SimdCap::AVX2,
        (int64_t)SimdCap::AVX512BW,
        (int64_t)SimdCap::AVX512VBMI2,
    }, {11, 12}, {32, 256, 4096}});

int main(int argc, char** argv) {
    init_simd_dispatch(nullptr);
    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
    return 0;
}
