#include "test_util.hpp"
#include "util/simd_dispatch.hpp"

#include <cstdio>
#include <cstdlib>
#include <string>

using namespace ikafssn;

static void test_init_no_logger() {
    reset_simd_dispatch_for_testing();
    init_simd_dispatch(nullptr);
    SimdCap cap = current_simd_cap();
    // Must be one of the known enum values (VBMI tier removed).
    bool valid = false;
    for (SimdCap c : {SimdCap::Scalar, SimdCap::SSE42, SimdCap::AVX2,
                      SimdCap::AVX512BW, SimdCap::AVX512VBMI2,
                      SimdCap::NEON, SimdCap::SVE, SimdCap::SVE2,
                      SimdCap::SME, SimdCap::SME2}) {
        if (cap == c) { valid = true; break; }
    }
    CHECK(valid);
}

static void test_simd_cap_name_complete() {
    // All known enum values must produce a name that is not "unknown".
    for (SimdCap c : {SimdCap::Scalar, SimdCap::SSE42, SimdCap::AVX2,
                      SimdCap::AVX512BW, SimdCap::AVX512VBMI2,
                      SimdCap::NEON, SimdCap::SVE, SimdCap::SVE2,
                      SimdCap::SME, SimdCap::SME2}) {
        CHECK(simd_cap_name(c) != "unknown");
    }
}

static void test_parse_force_simd_env_basic() {
    // Empty / null
    CHECK(parse_force_simd_env(nullptr).explicit_value == false);
    CHECK(parse_force_simd_env("").explicit_value == false);

    // Recognized tokens (case + separator insensitive). "scalar" is
    // recognised at the parser level even though init_simd_dispatch()
    // will exit(2) when it ends up as the active tier.
    CHECK(parse_force_simd_env("scalar").explicit_value);
    CHECK(parse_force_simd_env("scalar").cap == SimdCap::Scalar);

    CHECK(parse_force_simd_env("AVX2").explicit_value);
    CHECK(parse_force_simd_env("AVX2").cap == SimdCap::AVX2);

    CHECK(parse_force_simd_env("avx-2").explicit_value);
    CHECK(parse_force_simd_env("avx-2").cap == SimdCap::AVX2);

    CHECK(parse_force_simd_env("AVX_512_BW").explicit_value);
    CHECK(parse_force_simd_env("AVX_512_BW").cap == SimdCap::AVX512BW);

    CHECK(parse_force_simd_env("avx512vbmi2").explicit_value);
    CHECK(parse_force_simd_env("avx512vbmi2").cap == SimdCap::AVX512VBMI2);

    // standalone VBMI tier is removed; "avx512vbmi" silent-demotes
    // to AVX512BW so callers that hard-coded that string keep working.
    CHECK(parse_force_simd_env("avx512vbmi").explicit_value);
    CHECK(parse_force_simd_env("avx512vbmi").cap == SimdCap::AVX512BW);
    CHECK(parse_force_simd_env("AVX-512-VBMI").cap == SimdCap::AVX512BW);

    CHECK(parse_force_simd_env("sve2").explicit_value);
    CHECK(parse_force_simd_env("sve2").cap == SimdCap::SVE2);

    // Unrecognized
    CHECK(parse_force_simd_env("garbage").explicit_value == false);
    CHECK(parse_force_simd_env("avx9999").explicit_value == false);
}

static void test_parse_force_simd_env_auto() {
    // "auto" returns explicit_value=true with the per-arch maximum cap.
    auto r = parse_force_simd_env("auto");
    CHECK(r.explicit_value);
    // The cap should be at least the auto-detected value when on a supported
    // arch (so force_simd_cap() with this value clamps back to auto_cap).
    CHECK(static_cast<int>(r.cap) >= static_cast<int>(auto_detected_simd_cap()));
}

static void test_force_simd_cap_downgrade() {
    reset_simd_dispatch_for_testing();
    init_simd_dispatch(nullptr);
    SimdCap auto_cap = auto_detected_simd_cap();

    // Downgrade to scalar always succeeds.
    SimdCap applied = force_simd_cap(SimdCap::Scalar);
    CHECK(applied == SimdCap::Scalar);
    CHECK(current_simd_cap() == SimdCap::Scalar);

    // Re-apply auto (using force with the auto-detected value itself).
    applied = force_simd_cap(auto_cap);
    CHECK(applied == auto_cap);

    // Over-request collapses to scalar (use a cap deliberately higher than
    // any plausible auto-detected value on x86 by requesting an aarch64 tier).
#if defined(__x86_64__) || defined(__i386__)
    applied = force_simd_cap(SimdCap::SVE2);
    CHECK(applied == SimdCap::Scalar);
    CHECK(current_simd_cap() == SimdCap::Scalar);
#elif defined(__aarch64__)
    applied = force_simd_cap(SimdCap::AVX512VBMI2);
    CHECK(applied == SimdCap::Scalar);
    CHECK(current_simd_cap() == SimdCap::Scalar);
#endif
}

static void test_reset_and_re_init() {
    // Reset, force the lowest supported tier (SSE4.2 on x86, NEON on
    // arm), and verify that the env was honored.  "scalar" cannot be
    // used here because the dispatcher rejects Scalar at startup with
    // exit(2), which would terminate the unit binary.
#if defined(__x86_64__) || defined(__i386__)
    const char* low_tier = "sse42";
    SimdCap     low_cap  = SimdCap::SSE42;
#elif defined(__aarch64__)
    const char* low_tier = "neon";
    SimdCap     low_cap  = SimdCap::NEON;
#else
    const char* low_tier = nullptr;
    SimdCap     low_cap  = SimdCap::Scalar;
#endif
    if (low_tier) {
        reset_simd_dispatch_for_testing();
        setenv("IKAFSSN_FORCE_SIMD", low_tier, 1);
        init_simd_dispatch(nullptr);
        CHECK(current_simd_cap() == low_cap);

        // Detected (auto) value should still reflect the actual CPU.
        SimdCap auto_cap = auto_detected_simd_cap();
        bool valid = false;
        for (SimdCap c : {SimdCap::Scalar, SimdCap::SSE42, SimdCap::AVX2,
                          SimdCap::AVX512BW, SimdCap::AVX512VBMI2,
                          SimdCap::NEON, SimdCap::SVE, SimdCap::SVE2}) {
            if (auto_cap == c) { valid = true; break; }
        }
        CHECK(valid);

        unsetenv("IKAFSSN_FORCE_SIMD");
    }

    // Reset again and confirm default detection works without env override.
    reset_simd_dispatch_for_testing();
    init_simd_dispatch(nullptr);
}

static void test_test_require_tier_skip() {
    // When IKAFSSN_TEST_REQUIRE_TIER requests a tier higher than the actual CPU
    // supports, check_required_tier_or_skip() must call exit(77). We exercise
    // the "no skip needed" code path here (env unset → no-op).
    unsetenv("IKAFSSN_TEST_REQUIRE_TIER");
    check_required_tier_or_skip();  // must not exit

    // Set a tier guaranteed to be supported (Scalar) and verify no-op.
    setenv("IKAFSSN_TEST_REQUIRE_TIER", "scalar", 1);
    check_required_tier_or_skip();
    unsetenv("IKAFSSN_TEST_REQUIRE_TIER");
}

int main() {
    // Make sure this binary itself respects the SKIP protocol, since it is
    // registered with add_ikafssn_simd_test().
    init_simd_dispatch(nullptr);
    check_required_tier_or_skip();

    test_init_no_logger();
    test_simd_cap_name_complete();
    test_parse_force_simd_env_basic();
    test_parse_force_simd_env_auto();
    test_force_simd_cap_downgrade();
    test_reset_and_re_init();
    test_test_require_tier_skip();

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
