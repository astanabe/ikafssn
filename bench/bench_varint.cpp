// Phase 3a — varint batch decode SIMD micro-benchmark.
//
// Generates a stream of LEB128 varints with a typical posting-delta width
// distribution and benchmarks varint_decode_batch() per SIMD tier vs the
// scalar varint_decode() loop.
//
// Run: ./bench_varint --benchmark_format=console

#include "core/varint.hpp"
#include "index/varint_simd.hpp"
#include "util/simd_dispatch.hpp"
#include "bench_common.hpp"

#include <benchmark/benchmark.h>

#include <cstdint>
#include <random>
#include <vector>

using namespace ikafssn;

static std::vector<std::uint8_t> gen_varint_stream(std::size_t n_values,
                                                   std::mt19937& rng) {
    std::vector<std::uint8_t> buf;
    buf.reserve(n_values * 5);
    std::uint8_t tmp[5];
    for (std::size_t i = 0; i < n_values; i++) {
        std::uint32_t r = rng();
        std::uint32_t v;
        switch (r & 0x3) {
            case 0: v = r & 0x7Fu; break;       // 1-byte
            case 1: v = r & 0x3FFFu; break;     // 2-byte
            case 2: v = r & 0x1FFFFFu; break;   // 3-byte
            default: v = r; break;               // up to 5-byte
        }
        std::size_t n = varint_encode(v, tmp);
        for (std::size_t j = 0; j < n; j++) buf.push_back(tmp[j]);
    }
    return buf;
}

// Scalar baseline: call varint_decode in a loop.
static void BM_VarintDecodeScalarLoop(benchmark::State& state) {
    std::size_t n_values = static_cast<std::size_t>(state.range(0));
    std::mt19937 rng(42);
    auto stream = gen_varint_stream(n_values, rng);
    std::vector<std::uint32_t> out(n_values);

    for (auto _ : state) {
        const std::uint8_t* p = stream.data();
        const std::uint8_t* e = stream.data() + stream.size();
        std::size_t i = 0;
        while (p < e && i < n_values) {
            std::uint32_t v;
            p += varint_decode(p, v);
            out[i++] = v;
        }
        benchmark::DoNotOptimize(out.data());
    }
    state.SetItemsProcessed(int64_t(state.iterations()) * n_values);
    state.SetBytesProcessed(int64_t(state.iterations()) * stream.size());
}
BENCHMARK(BM_VarintDecodeScalarLoop)->Arg(1 << 10)->Arg(1 << 14)->Arg(1 << 18);

// Batch-decode at a specified tier.
static void BM_VarintDecodeBatch(benchmark::State& state) {
    SimdCap tier = static_cast<SimdCap>(state.range(0));
    int batch = static_cast<int>(state.range(1));
    std::size_t n_values = static_cast<std::size_t>(state.range(2));
    bench::apply_force_tier(state, tier);
    bench::label_tier(state, tier);

    std::mt19937 rng(42);
    auto stream = gen_varint_stream(n_values, rng);
    std::vector<std::uint32_t> out(n_values);

    for (auto _ : state) {
        const std::uint8_t* p = stream.data();
        const std::uint8_t* e = stream.data() + stream.size();
        std::size_t produced = 0;
        while (p < e && produced < n_values) {
            std::size_t consumed = 0;
            int n = varint_decode_batch(p, e, out.data() + produced,
                                        &consumed, batch);
            if (n == 0) break;
            p += consumed;
            produced += n;
        }
        benchmark::DoNotOptimize(out.data());
    }
    state.SetItemsProcessed(int64_t(state.iterations()) * n_values);
    state.SetBytesProcessed(int64_t(state.iterations()) * stream.size());
}
BENCHMARK(BM_VarintDecodeBatch)
    ->ArgsProduct({{
        (int64_t)SimdCap::Scalar,
        (int64_t)SimdCap::SSE42,
        (int64_t)SimdCap::AVX2,
        (int64_t)SimdCap::AVX512BW,
        (int64_t)SimdCap::AVX512VBMI,
        (int64_t)SimdCap::AVX512VBMI2,
    }, {16}, {1 << 10, 1 << 14, 1 << 18}});

int main(int argc, char** argv) {
    init_simd_dispatch(nullptr);
    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
    return 0;
}
