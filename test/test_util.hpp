#pragma once

#include <cstdio>
#include <cstdlib>

static int g_fail_count = 0;
static int g_pass_count = 0;

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
