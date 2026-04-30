#pragma once

#include <cstdio>
#include <cstdlib>
#include <string>

static int g_fail_count = 0;
static int g_pass_count = 0;

// Build a per-process temp directory path. When the test is invoked as a
// SIMD variant (`add_ikafssn_simd_test` sets IKAFSSN_FORCE_SIMD), the tier
// name is appended so that the base test and every variant have disjoint
// directories and can run truly in parallel under `ctest -j`.
inline std::string test_tmpdir(const char* base) {
    const char* simd = std::getenv("IKAFSSN_FORCE_SIMD");
    if (simd && *simd) {
        return std::string(base) + "_" + simd;
    }
    return std::string(base);
}

#define CHECK(cond) \
    do { \
        if (!(cond)) { \
            std::fprintf(stderr, "FAIL: %s:%d: %s\n", __FILE__, __LINE__, #cond); \
            g_fail_count++; \
        } else { \
            g_pass_count++; \
        } \
    } while (0)

#define CHECK_EQ(a, b) \
    do { \
        auto _a = (a); \
        auto _b = (b); \
        if (_a != _b) { \
            std::fprintf(stderr, "FAIL: %s:%d: %s == %s  (got %llu vs %llu)\n", \
                         __FILE__, __LINE__, #a, #b, \
                         (unsigned long long)_a, (unsigned long long)_b); \
            g_fail_count++; \
        } else { \
            g_pass_count++; \
        } \
    } while (0)

#define TEST_SUMMARY() \
    do { \
        std::fprintf(stderr, "%d passed, %d failed\n", g_pass_count, g_fail_count); \
    } while (0)
