// Phase 5a — prototype micro-benchmark for SIMD-FastPFOR* + Simple-8b decode.
//
// Generates a stream of synthetic uint32 deltas matching the typical posting
// distribution (mostly 1- or 2-byte after delta encoding, with a long tail)
// and benchmarks:
//
//   1. Phase 3a varint scalar-loop decode (baseline)
//   2. Phase 3a varint_decode_batch SIMD decode (current production)
//   3. SIMD-FastPFOR* + VariableByte tail decode (FastPFOR 0.4.0, codec
//      name "simdfastpfor256") — proxy for proposed v4 .kix codec
//
// Goal: produce a speed_ratio numeric for the Phase 5a go/no-go gate
// (size_ratio <= 0.90 OR speed_ratio <= 0.90).
//
// Run: ./bench_pfd_proto --benchmark_format=console

#include "core/varint.hpp"
#include "index/varint_simd.hpp"
#include "util/simd_dispatch.hpp"
#include "bench_common.hpp"

#include <benchmark/benchmark.h>

// FastPFOR is added via FetchContent in the parent CMakeLists.txt.
#include "codecfactory.h"

#include <cstdint>
#include <memory>
#include <random>
#include <vector>

using namespace ikafssn;

namespace {

// Same delta distribution as bench_varint.cpp so numbers are comparable.
std::vector<std::uint32_t> gen_delta_stream(std::size_t n_values,
                                            std::mt19937& rng) {
    std::vector<std::uint32_t> out;
    out.reserve(n_values);
    for (std::size_t i = 0; i < n_values; i++) {
        std::uint32_t r = rng();
        std::uint32_t v;
        switch (r & 0x3) {
            case 0: v = r & 0x7Fu; break;       // 1-byte varint
            case 1: v = r & 0x3FFFu; break;     // 2-byte varint
            case 2: v = r & 0x1FFFFFu; break;   // 3-byte varint
            default: v = r; break;               // up to 5-byte varint
        }
        out.push_back(v);
    }
    return out;
}

std::vector<std::uint8_t> encode_varints(const std::vector<std::uint32_t>& deltas) {
    std::vector<std::uint8_t> buf;
    buf.reserve(deltas.size() * 5);
    std::uint8_t tmp[5];
    for (auto d : deltas) {
        std::size_t n = varint_encode(d, tmp);
        for (std::size_t j = 0; j < n; j++) buf.push_back(tmp[j]);
    }
    return buf;
}

// Encode using FastPFOR's CompositeCodec<SIMDFastPFor<8>, VariableByte> codec.
// The codec writes to a uint32_t buffer (it's a 4-byte-granular layout).
std::vector<std::uint32_t> encode_pfor(FastPForLib::IntegerCODEC& codec,
                                       const std::vector<std::uint32_t>& deltas) {
    std::vector<std::uint32_t> out(deltas.size() + 1024);
    std::size_t nvalue = out.size();
    codec.encodeArray(deltas.data(), deltas.size(), out.data(), nvalue);
    out.resize(nvalue);
    return out;
}

} // namespace

// === Baseline: scalar varint loop ===
static void BM_PfdProto_VarintScalar(benchmark::State& state) {
    std::size_t n_values = static_cast<std::size_t>(state.range(0));
    std::mt19937 rng(42);
    auto deltas = gen_delta_stream(n_values, rng);
    auto stream = encode_varints(deltas);
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
BENCHMARK(BM_PfdProto_VarintScalar)->Arg(1 << 14)->Arg(1 << 18);

// === Phase 3a SIMD varint batch ===
static void BM_PfdProto_VarintSimd(benchmark::State& state) {
    SimdCap tier = static_cast<SimdCap>(state.range(0));
    std::size_t n_values = static_cast<std::size_t>(state.range(1));
    bench::apply_force_tier(state, tier);
    bench::label_tier(state, tier);

    std::mt19937 rng(42);
    auto deltas = gen_delta_stream(n_values, rng);
    auto stream = encode_varints(deltas);
    std::vector<std::uint32_t> out(n_values);

    constexpr int batch = 16;
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
BENCHMARK(BM_PfdProto_VarintSimd)
    ->ArgsProduct({{
        (int64_t)SimdCap::AVX2,
        (int64_t)SimdCap::AVX512BW,
        (int64_t)SimdCap::AVX512VBMI,
        (int64_t)SimdCap::AVX512VBMI2,
    }, {1 << 14, 1 << 18}});

// === SIMD-FastPFOR* + VariableByte tail ("simdfastpfor256") ===
//
// This is the closest factory codec to the plan's "SIMD-FastPFOR* +
// Simple-8b" layout.  It uses 256-elem blocks rather than 128 (Phase 5b
// will use 128), but the per-int decode cost is the relevant signal.
static void BM_PfdProto_FastPforDecode(benchmark::State& state) {
    std::size_t n_values = static_cast<std::size_t>(state.range(0));
    std::mt19937 rng(42);
    auto deltas = gen_delta_stream(n_values, rng);

    static FastPForLib::CODECFactory factory;
    auto codec = factory.getFromName("simdfastpfor256");
    if (!codec) {
        state.SkipWithError("simdfastpfor256 codec not found");
        return;
    }
    auto encoded = encode_pfor(*codec, deltas);
    state.counters["encoded_bytes"] = encoded.size() * 4.0;
    state.counters["raw_bytes"] = double(deltas.size() * 4);

    std::vector<std::uint32_t> out(n_values + 1024);
    for (auto _ : state) {
        std::size_t nvalue = out.size();
        codec->decodeArray(encoded.data(), encoded.size(), out.data(), nvalue);
        benchmark::DoNotOptimize(out.data());
    }
    state.SetItemsProcessed(int64_t(state.iterations()) * n_values);
    state.SetBytesProcessed(int64_t(state.iterations()) * encoded.size() * 4);
}
BENCHMARK(BM_PfdProto_FastPforDecode)->Arg(1 << 14)->Arg(1 << 18);

// === Simple8b alone (tail codec proxy) ===
static void BM_PfdProto_Simple8bDecode(benchmark::State& state) {
    std::size_t n_values = static_cast<std::size_t>(state.range(0));
    std::mt19937 rng(42);
    auto deltas = gen_delta_stream(n_values, rng);

    static FastPForLib::CODECFactory factory;
    auto codec = factory.getFromName("simple8b");
    if (!codec) {
        state.SkipWithError("simple8b codec not found");
        return;
    }
    auto encoded = encode_pfor(*codec, deltas);
    state.counters["encoded_bytes"] = encoded.size() * 4.0;
    state.counters["raw_bytes"] = double(deltas.size() * 4);

    std::vector<std::uint32_t> out(n_values + 1024);
    for (auto _ : state) {
        std::size_t nvalue = out.size();
        codec->decodeArray(encoded.data(), encoded.size(), out.data(), nvalue);
        benchmark::DoNotOptimize(out.data());
    }
    state.SetItemsProcessed(int64_t(state.iterations()) * n_values);
    state.SetBytesProcessed(int64_t(state.iterations()) * encoded.size() * 4);
}
BENCHMARK(BM_PfdProto_Simple8bDecode)->Arg(1 << 14)->Arg(1 << 18);

int main(int argc, char** argv) {
    init_simd_dispatch(nullptr);
    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
    return 0;
}
