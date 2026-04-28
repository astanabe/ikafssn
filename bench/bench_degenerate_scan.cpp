// Phase 4a — degenerate-base scan SIMD micro-benchmark.
//
// Measures has_degenerate_base() at each available SIMD tier, plus a scalar
// baseline that mirrors the pre-Phase-4a inline body. The "clean" payload
// (pure ACGT) exercises the chunk loop to completion; "head_hit" / "tail_hit"
// exercise the early-exit branches.

#include "core/degenerate_scan_simd.hpp"
#include "core/kmer_encoding.hpp"
#include "util/simd_dispatch.hpp"
#include "bench_common.hpp"

#include <benchmark/benchmark.h>

#include <cstddef>
#include <cstdint>
#include <random>
#include <string>
#include <vector>

using namespace ikafssn;

enum Payload {
    kPayloadClean    = 0,
    kPayloadHeadHit  = 1,
    kPayloadTailHit  = 2,
};

static std::string make_payload(std::size_t n, Payload kind, std::mt19937& rng) {
    std::string s(n, 'A');
    const char ab[4] = {'A', 'C', 'G', 'T'};
    for (std::size_t i = 0; i < n; i++) s[i] = ab[rng() & 3];
    if (kind == kPayloadHeadHit && n > 0) s[0] = 'N';
    if (kind == kPayloadTailHit && n > 0) s[n - 1] = 'R';
    return s;
}

// Scalar baseline: per-byte LUT lookup (pre-Phase-4a behavior).
static void BM_DegenerateScanScalarLoop(benchmark::State& state) {
    std::size_t n = static_cast<std::size_t>(state.range(0));
    Payload kind  = static_cast<Payload>(state.range(1));
    std::mt19937 rng(42);
    std::string s = make_payload(n, kind, rng);
    const bool* tbl = degenerate_base_table();
    for (auto _ : state) {
        bool hit = false;
        for (char c : s) {
            if (tbl[static_cast<std::uint8_t>(c)]) { hit = true; break; }
        }
        benchmark::DoNotOptimize(hit);
    }
    state.SetItemsProcessed(int64_t(state.iterations()) * n);
    state.SetBytesProcessed(int64_t(state.iterations()) * n);
}
BENCHMARK(BM_DegenerateScanScalarLoop)
    ->ArgsProduct({{32, 256, 1024, 4096, 16384},
                   {kPayloadClean, kPayloadHeadHit, kPayloadTailHit}});

static void BM_DegenerateScan(benchmark::State& state) {
    SimdCap tier  = static_cast<SimdCap>(state.range(0));
    std::size_t n = static_cast<std::size_t>(state.range(1));
    Payload kind  = static_cast<Payload>(state.range(2));
    bench::apply_force_tier(state, tier);
    bench::label_tier(state, tier);

    std::mt19937 rng(42);
    std::string s = make_payload(n, kind, rng);
    for (auto _ : state) {
        bool hit = has_degenerate_base(s.data(), s.size());
        benchmark::DoNotOptimize(hit);
    }
    state.SetItemsProcessed(int64_t(state.iterations()) * n);
    state.SetBytesProcessed(int64_t(state.iterations()) * n);
}
BENCHMARK(BM_DegenerateScan)
    ->ArgsProduct({{
        (int64_t)SimdCap::SSE42,
        (int64_t)SimdCap::AVX2,
        (int64_t)SimdCap::AVX512BW,
        (int64_t)SimdCap::AVX512VBMI2,
    }, {32, 256, 1024, 4096, 16384},
       {kPayloadClean, kPayloadHeadHit, kPayloadTailHit}});

int main(int argc, char** argv) {
    init_simd_dispatch(nullptr);
    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
    return 0;
}
