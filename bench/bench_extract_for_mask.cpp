// extract_for_mask_batch SIMD micro-benchmark.
//
// Measures the batch entry point at each available SIMD tier vs a
// scalar baseline that mirrors the per-element extract_for_mask()
// call site.

#include "core/spaced_seed.hpp"
#include "core/spaced_seed_simd.hpp"
#include "util/simd_dispatch.hpp"
#include "bench_common.hpp"

#include <benchmark/benchmark.h>

#include <cstdint>
#include <random>
#include <vector>

using namespace ikafssn;

namespace {

// 6 representative masks across all (k, t, type) families.
struct MaskInfo {
    std::uint32_t mask;
    bool          is_u32;   // false = uint16_t, true = uint32_t
    const char*   label;
};

constexpr MaskInfo kMasks[] = {
    {MASK_K8_T13_CODING,   false, "K8_T13_COD"},
    {MASK_K8_T18_CODING,   false, "K8_T18_COD"},
    {MASK_K9_T18_OPTIMAL,  true,  "K9_T18_OPT"},
    {MASK_K11_T16_CODING,  true,  "K11_T16_COD"},
    {MASK_K11_T21_CODING,  true,  "K11_T21_COD"},
    {MASK_K12_T18_OPTIMAL, true,  "K12_T18_OPT"},
};

template <typename KmerInt>
std::vector<std::uint64_t> make_input(std::size_t n, std::mt19937& rng) {
    std::vector<std::uint64_t> a(n);
    for (auto& v : a) {
        v = (std::uint64_t)rng() ^ ((std::uint64_t)rng() << 32);
    }
    return a;
}

} // namespace

// Scalar baseline: per-element extract_for_mask() in a loop.
template <typename KmerInt>
static void BM_ExtractForMaskScalarLoop(benchmark::State& state) {
    int mask_idx = static_cast<int>(state.range(0));
    std::size_t n = static_cast<std::size_t>(state.range(1));
    std::uint32_t mask = kMasks[mask_idx].mask;
    state.SetLabel(kMasks[mask_idx].label);
    std::mt19937 rng(42);
    auto a = make_input<KmerInt>(n, rng);
    std::vector<KmerInt> out(n);
    for (auto _ : state) {
        for (std::size_t i = 0; i < n; ++i) {
            out[i] = extract_for_mask<KmerInt>(a[i], mask);
        }
        benchmark::DoNotOptimize(out.data());
    }
    state.SetItemsProcessed(int64_t(state.iterations()) * n);
}
// uint16_t masks (indices 0, 1)
BENCHMARK(BM_ExtractForMaskScalarLoop<std::uint16_t>)
    ->ArgsProduct({{0, 1}, {16, 64, 256, 1024, 4096}});
// uint32_t masks (indices 2..5)
BENCHMARK(BM_ExtractForMaskScalarLoop<std::uint32_t>)
    ->ArgsProduct({{2, 3, 4, 5}, {16, 64, 256, 1024, 4096}});

template <typename KmerInt>
static void BM_ExtractForMaskBatch(benchmark::State& state) {
    SimdCap tier  = static_cast<SimdCap>(state.range(0));
    int mask_idx  = static_cast<int>(state.range(1));
    std::size_t n = static_cast<std::size_t>(state.range(2));
    std::uint32_t mask = kMasks[mask_idx].mask;
    bench::apply_force_tier(state, tier);
    state.SetLabel(std::string(simd_cap_name(tier)) + "/" + kMasks[mask_idx].label);

    std::mt19937 rng(42);
    auto a = make_input<KmerInt>(n, rng);
    std::vector<KmerInt> out(n);
    for (auto _ : state) {
        extract_for_mask_batch<KmerInt>(a.data(), n, mask, out.data());
        benchmark::DoNotOptimize(out.data());
    }
    state.SetItemsProcessed(int64_t(state.iterations()) * n);
}
BENCHMARK(BM_ExtractForMaskBatch<std::uint16_t>)
    ->ArgsProduct({{
        (int64_t)SimdCap::SSE42,
        (int64_t)SimdCap::AVX2,
        (int64_t)SimdCap::AVX512BW,
        (int64_t)SimdCap::AVX512VBMI2,
    }, {0, 1}, {16, 64, 256, 1024, 4096}});
BENCHMARK(BM_ExtractForMaskBatch<std::uint32_t>)
    ->ArgsProduct({{
        (int64_t)SimdCap::SSE42,
        (int64_t)SimdCap::AVX2,
        (int64_t)SimdCap::AVX512BW,
        (int64_t)SimdCap::AVX512VBMI2,
    }, {2, 3, 4, 5}, {16, 64, 256, 1024, 4096}});

int main(int argc, char** argv) {
    init_simd_dispatch(nullptr);
    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
    return 0;
}
