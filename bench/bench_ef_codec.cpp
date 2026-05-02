// Phase 7b — Elias-Fano dictionary codec micro-benchmark.
//
// Measures encode_dictionary_ef throughput and EFDictionary::access
// per-call cost across the runtime SIMD ladder (sse42 / avx2 /
// avx512bw / avx512vbmi2 on x86; neon on aarch64).  The hot-path
// access() figure is the eliasfano.md §7 driver — at AVX2-and-above
// tiers the inner-word select1 uses BMI2 PDEP+TZCNT; below that it
// falls back to __builtin_popcountll + bit-stripping __builtin_ctzll.

#include "index/ef_codec.hpp"
#include "index/ef_format.hpp"
#include "util/simd_dispatch.hpp"
#include "bench_common.hpp"

#include <benchmark/benchmark.h>

#include <algorithm>
#include <cstdint>
#include <random>
#include <vector>

using namespace ikafssn;

namespace {

// Build a synthetic monotonically-non-decreasing offsets array of
// length D where the average gap is U/D (i.e. l = log2(U/D)).  Mirrors
// the .kix / .kpx dictionary shape.
std::vector<std::uint64_t> build_offsets(std::size_t D, std::uint64_t U,
                                          std::uint32_t seed) {
    std::vector<std::uint64_t> offsets;
    offsets.reserve(D);
    std::mt19937 rng(seed);
    if (D == 0) return offsets;
    std::uint64_t avg_step = std::max<std::uint64_t>(1, U / D);
    std::uint64_t off = 0;
    for (std::size_t i = 0; i < D; ++i) {
        offsets.push_back(off);
        // 25 % chance of a "run" (offset stays the same — empty
        // posting list); otherwise advance by a random fraction of
        // the average step to keep total under U.
        if ((rng() & 3) != 0) {
            off += rng() % (2 * avg_step + 1);
        }
    }
    return offsets;
}

} // anonymous namespace

// ===== EF encode throughput =====

static void BM_EFEncode(benchmark::State& state) {
    SimdCap tier = static_cast<SimdCap>(state.range(0));
    std::size_t D = static_cast<std::size_t>(state.range(1));
    bench::apply_force_tier(state, tier);
    bench::label_tier(state, tier);

    const std::uint64_t U = static_cast<std::uint64_t>(D) * 16ULL;
    auto offsets = build_offsets(D, U, 0xC0FFEEu);
    std::vector<std::uint8_t> blob;
    blob.reserve(64 * 1024);

    for (auto _ : state) {
        blob.clear();
        ef::encode_dictionary_ef(offsets.data(), D, U, blob);
        benchmark::DoNotOptimize(blob.data());
    }
    state.SetItemsProcessed(int64_t(state.iterations()) * D);
    state.SetBytesProcessed(int64_t(state.iterations()) * D * sizeof(std::uint64_t));
}
BENCHMARK(BM_EFEncode)
    ->ArgsProduct({{
        (int64_t)SimdCap::SSE42,
        (int64_t)SimdCap::AVX2,
        (int64_t)SimdCap::AVX512BW,
        (int64_t)SimdCap::AVX512VBMI2,
    }, {1024, 16384, 262144}});

// ===== EF access (random) — Stage 1 hot path =====

static void BM_EFAccessRandom(benchmark::State& state) {
    SimdCap tier = static_cast<SimdCap>(state.range(0));
    std::size_t D = static_cast<std::size_t>(state.range(1));
    bench::apply_force_tier(state, tier);
    bench::label_tier(state, tier);

    const std::uint64_t U = static_cast<std::uint64_t>(D) * 16ULL;
    auto offsets = build_offsets(D, U, 0xBEEFu);
    std::vector<std::uint8_t> blob;
    ef::encode_dictionary_ef(offsets.data(), D, U, blob);

    ef::EFDictionary dict;
    if (!dict.open(blob.data(), blob.size())) {
        state.SkipWithError("EFDictionary::open failed");
        return;
    }

    std::mt19937 rng(0xCAFEu);
    std::vector<std::uint32_t> indices(4096);
    for (auto& v : indices) v = static_cast<std::uint32_t>(rng() % D);

    std::uint64_t sink = 0;
    for (auto _ : state) {
        for (std::uint32_t i : indices) {
            sink += dict.access(i);
        }
        benchmark::DoNotOptimize(sink);
    }
    state.SetItemsProcessed(int64_t(state.iterations()) * indices.size());
}
BENCHMARK(BM_EFAccessRandom)
    ->ArgsProduct({{
        (int64_t)SimdCap::SSE42,
        (int64_t)SimdCap::AVX2,
        (int64_t)SimdCap::AVX512BW,
        (int64_t)SimdCap::AVX512VBMI2,
    }, {1024, 16384, 262144, 4194304}});

// ===== EF access (sequential) — bulk_count_postings shape =====

static void BM_EFAccessSequential(benchmark::State& state) {
    SimdCap tier = static_cast<SimdCap>(state.range(0));
    std::size_t D = static_cast<std::size_t>(state.range(1));
    bench::apply_force_tier(state, tier);
    bench::label_tier(state, tier);

    const std::uint64_t U = static_cast<std::uint64_t>(D) * 16ULL;
    auto offsets = build_offsets(D, U, 0xFEEDu);
    std::vector<std::uint8_t> blob;
    ef::encode_dictionary_ef(offsets.data(), D, U, blob);

    ef::EFDictionary dict;
    if (!dict.open(blob.data(), blob.size())) {
        state.SkipWithError("EFDictionary::open failed");
        return;
    }

    std::uint64_t sink = 0;
    for (auto _ : state) {
        std::uint64_t prev = dict.access(0);
        for (std::size_t i = 1; i < D; ++i) {
            std::uint64_t cur = dict.access(static_cast<std::uint32_t>(i));
            sink += cur - prev;
            prev = cur;
        }
        benchmark::DoNotOptimize(sink);
    }
    state.SetItemsProcessed(int64_t(state.iterations()) * D);
}
BENCHMARK(BM_EFAccessSequential)
    ->ArgsProduct({{
        (int64_t)SimdCap::SSE42,
        (int64_t)SimdCap::AVX2,
        (int64_t)SimdCap::AVX512BW,
        (int64_t)SimdCap::AVX512VBMI2,
    }, {1024, 16384, 262144}});

int main(int argc, char** argv) {
    init_simd_dispatch(nullptr);
    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
    return 0;
}
