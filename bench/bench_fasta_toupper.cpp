#include "io/text_simd.hpp"
#include "util/simd_dispatch.hpp"
#include "bench/bench_common.hpp"

#include <benchmark/benchmark.h>

#include <cstdint>
#include <random>
#include <vector>

using namespace ikafssn;

static void BM_ToupperAscii(benchmark::State& state) {
    SimdCap     tier = static_cast<SimdCap>(state.range(0));
    std::size_t n    = static_cast<std::size_t>(state.range(1));
    bench::apply_force_tier(state, tier);
    bench::label_tier(state, tier);

    std::vector<std::uint8_t> buf(n);
    std::mt19937 rng(42);
    for (auto& b : buf) b = static_cast<std::uint8_t>('a' + (rng() % 26));

    for (auto _ : state) {
        toupper_inplace_ascii(buf.data(), buf.size());
        // Re-lowercase a sparse subset so subsequent iterations still have
        // work to do (otherwise the loop becomes a no-op after iter 1).
        state.PauseTiming();
        for (std::size_t i = 0; i < n; i += 4096) buf[i] |= 0x20;
        state.ResumeTiming();
    }
    state.SetBytesProcessed(static_cast<std::int64_t>(state.iterations()) *
                            static_cast<std::int64_t>(n));
}

static void RegisterToupperArgs(benchmark::internal::Benchmark* b) {
    static constexpr int kTiers[] = {
        static_cast<int>(SimdCap::Scalar),
#if defined(__x86_64__) || defined(__i386__)
        static_cast<int>(SimdCap::SSE42),
        static_cast<int>(SimdCap::AVX2),
        static_cast<int>(SimdCap::AVX512BW),
        static_cast<int>(SimdCap::AVX512VBMI),
        static_cast<int>(SimdCap::AVX512VBMI2),
#elif defined(__aarch64__)
        static_cast<int>(SimdCap::NEON),
        static_cast<int>(SimdCap::SVE),
        static_cast<int>(SimdCap::SVE2),
#endif
    };
    static constexpr std::int64_t kSizes[] = {
        64, 1024, 65536, 1 << 20, 16 << 20
    };
    for (int t : kTiers) {
        for (std::int64_t n : kSizes) b->Args({t, n});
    }
}

BENCHMARK(BM_ToupperAscii)->Apply(RegisterToupperArgs)->UseRealTime();

int main(int argc, char** argv) {
    init_simd_dispatch(nullptr);
    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
    return 0;
}
