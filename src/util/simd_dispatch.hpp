#pragma once

#include "util/logger.hpp"

#include <cstddef>
#include <cstdint>
#include <string_view>

namespace ikafssn {

// Runtime SIMD capability tier. Values are arch-segregated by range so that
// cross-arch comparisons are structurally no-ops:
//   x86_64 :  10..50
//   aarch64: 100..199
enum class SimdCap : int {
    // Scalar = 0 is the detection-time sentinel (no SIMD detected) and
    // the over-request collapse target for force_simd_cap().  SSE4.2
    // is the production x86_64 floor: init_simd_dispatch() rejects
    // auto_cap == Scalar with exit(2) unless the build was configured
    // with IKAFSSN_ENABLE_SIMD=0 (build_disabled path).
    Scalar       =   0,
    // x86_64
    SSE42        =  10,
    AVX2         =  20,  // BMI1/BMI2 implied; slow_bmi2 is tracked separately
    AVX512BW     =  30,  // standalone VBMI is silently demoted to this tier
    AVX512VBMI2  =  50,
    // aarch64
    NEON         = 100,
    SVE          = 110,
    SVE2         = 120,
    SME          = 130,  // detected only; not used by current kernels
    SME2         = 140,  // detected only; not used by current kernels
};

struct SimdAuxFlags {
    bool slow_bmi2        = false;  // AMD Zen1/Zen2 (microcoded pdep/pext)
    bool has_avx512f_only = false;  // KNL/KNM (demoted to AVX2 tier; informational)
};

// Advisory minimum input size (bytes) for SIMD kernels at each tier.
inline constexpr std::size_t kSimdMinBytes_SSE42       =  16;
inline constexpr std::size_t kSimdMinBytes_AVX2        =  64;
inline constexpr std::size_t kSimdMinBytes_AVX512BW    = 128;
inline constexpr std::size_t kSimdMinBytes_AVX512VBMI2 = 512;
inline constexpr std::size_t kSimdMinBytes_NEON        =  16;
inline constexpr std::size_t kSimdMinBytes_SVE         =  64;
inline constexpr std::size_t kSimdMinBytes_SVE2        =  64;

// Initialize the runtime SIMD dispatch state. Idempotent: protected by
// std::call_once. Safe to call from multiple binaries / threads. The logger
// pointer may be nullptr for silent operation.
void init_simd_dispatch(Logger* logger = nullptr) noexcept;

// Currently active capability (may be downgraded by IKAFSSN_FORCE_SIMD or by
// IKAFSSN_ENABLE_SIMD=0 build-time disable).
[[nodiscard]] SimdCap         current_simd_cap()       noexcept;

// CPU-detected capability before any IKAFSSN_FORCE_SIMD downgrade.
[[nodiscard]] SimdCap         auto_detected_simd_cap() noexcept;

[[nodiscard]] std::string_view simd_cap_name(SimdCap)  noexcept;

// Bench/fixture entrypoint: requested tier must be at most the auto-detected
// value. Over-requests collapse to Scalar to prevent SIGILL on unsupported CPUs.
// Returns the actually applied capability.
SimdCap force_simd_cap(SimdCap requested) noexcept;

struct ParseForceResult {
    SimdCap          cap;
    bool             explicit_value;  // true = recognized token (incl. "auto"/"scalar")
    std::string_view token;            // diagnostic echo
};
[[nodiscard]] ParseForceResult parse_force_simd_env(const char* env_value) noexcept;

// CTest variant skip helper: if env var IKAFSSN_TEST_REQUIRE_TIER names a tier
// that the actual CPU (auto-detected, not the possibly-downgraded current tier)
// does not support, exit(77) so CTest records a SKIP. Cross-arch tier requests
// are also skipped (auto-detect range guarantees the comparison fails).
// Must be called after init_simd_dispatch().
void check_required_tier_or_skip() noexcept;

#ifdef IKAFSSN_TESTING
// Re-arm the call_once for unit-test reset. Not safe in production code.
void reset_simd_dispatch_for_testing() noexcept;
#endif

} // namespace ikafssn
