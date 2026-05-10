#pragma once

#include "util/simd_dispatch.hpp"

#include <benchmark/benchmark.h>

#include <string>

namespace ikafssn::bench {

// Force a SIMD tier inside a benchmark fixture. force_simd_cap() is
// downgrade-only, so requesting a tier above the auto-detected cap
// collapses to Scalar — surface that as a SkipWithError instead of
// silently producing a misleading datapoint.
inline void apply_force_tier(benchmark::State& state, SimdCap requested) {
    SimdCap got = force_simd_cap(requested);
    if (got != requested) {
        state.SkipWithError(
            "force_simd_cap downgraded: CPU does not support requested tier");
    }
}

inline void label_tier(benchmark::State& state, SimdCap tier) {
    state.SetLabel(std::string(simd_cap_name(tier)));
}

} // namespace ikafssn::bench
