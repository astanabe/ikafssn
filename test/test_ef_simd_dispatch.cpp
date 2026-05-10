#include "test_util.hpp"
#include "index/ef_codec.hpp"
#include "util/simd_dispatch.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <string_view>

using namespace ikafssn;

namespace {

void test_active_tier_matches_force_simd() {
    const char* forced = std::getenv("IKAFSSN_FORCE_SIMD");
    if (!forced || *forced == '\0') {
        // Without IKAFSSN_FORCE_SIMD, just check the dispatcher returns
        // a recognised tier name.
        std::string_view name = ef::active_tier_name();
        bool valid = (name == "sse42") || (name == "avx2")
                  || (name == "avx512bw") || (name == "avx512vbmi2")
                  || (name == "neon");
        CHECK(valid);
        return;
    }

    // silent demotions: avx512vbmi → avx512bw.
    std::string expected = forced;
    if (expected == "avx512vbmi" || expected == "AVX512VBMI") {
        expected = "avx512bw";
    }

    std::string_view active = ef::active_tier_name();
    if (active != expected) {
        std::fprintf(stderr,
            "FAIL: ef::active_tier_name() = %s, expected %s "
            "(IKAFSSN_FORCE_SIMD=%s)\n",
            std::string(active).c_str(), expected.c_str(), forced);
    }
    CHECK(active == expected);
}

} // anonymous namespace

int main() {
    init_simd_dispatch(nullptr);
    check_required_tier_or_skip();

    test_active_tier_matches_force_simd();

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
