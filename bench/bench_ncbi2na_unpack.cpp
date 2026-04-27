#include "core/ncbi2na_unpack.hpp"
#include "util/simd_dispatch.hpp"
#include "bench/bench_common.hpp"

#include <benchmark/benchmark.h>

#include <cstdint>
#include <random>
#include <vector>

using namespace ikafssn;

static void BM_Ncbi2naUnpack(benchmark::State& state) {
    SimdCap     tier    = static_cast<SimdCap>(state.range(0));
    std::size_t n_bases = static_cast<std::size_t>(state.range(1));
    bench::apply_force_tier(state, tier);
    bench::label_tier(state, tier);

    std::vector<char>          packed((n_bases + 3) / 4 + 16);  // headroom
    std::vector<std::uint8_t>  out(n_bases);
    std::mt19937 rng(42);
    for (auto& b : packed) b = static_cast<char>(rng() & 0xFF);

    for (auto _ : state) {
        for (std::size_t pos = 0; pos + 64 <= n_bases; pos += 64) {
            unpack_ncbi2na_chunk(packed.data(),
                                 static_cast<std::uint32_t>(pos),
                                 64,
                                 out.data() + pos);
        }
        benchmark::DoNotOptimize(out.data());
    }
    state.SetBytesProcessed(static_cast<std::int64_t>(state.iterations()) *
                            static_cast<std::int64_t>(n_bases));
}

static void RegisterUnpackArgs(benchmark::internal::Benchmark* b) {
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
        256, 4096, 65536, 1 << 20, 16 << 20
    };
    for (int t : kTiers) {
        for (std::int64_t n : kSizes) b->Args({t, n});
    }
}

BENCHMARK(BM_Ncbi2naUnpack)->Apply(RegisterUnpackArgs)->UseRealTime();

int main(int argc, char** argv) {
    init_simd_dispatch(nullptr);
    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
    return 0;
}
