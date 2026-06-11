#include "test_util.hpp"
#include "index/pfd_codec.hpp"

#include <cstdint>
#include <cstdio>
#include <random>
#include <vector>

using namespace ikafssn;

namespace {

// Encode a 2-bit kind for index `i` into `kind_map`.  Mirrors the
// set_kind_bits helper in pfd_codec_tier.cpp; replicated here so the
// test does not depend on internal symbols.
inline void set_kind(std::uint8_t* kind_map, std::uint32_t i, std::uint8_t kind) {
    kind_map[i >> 2] |= std::uint8_t((kind & 0x03) << ((i & 3) * 2));
}

void check_popcount(const std::vector<std::uint8_t>& kinds,
                    const char* label) {
    const std::uint32_t n = static_cast<std::uint32_t>(kinds.size());
    const std::size_t bytes = (std::size_t(n) * 2 + 7) / 8;
    std::vector<std::uint8_t> km(bytes, 0);
    std::uint32_t expected_partition = 0, expected_short1 = 0, expected_short2 = 0;
    for (std::uint32_t i = 0; i < n; ++i) {
        set_kind(km.data(), i, kinds[i]);
        switch (kinds[i]) {
            case 0: expected_short1++;    break;
            case 1: expected_short2++;    break;
            case 2: expected_partition++; break;
        }
    }

    std::uint32_t got_partition = 0xFFFFFFFFu;
    std::uint32_t got_short1    = 0xFFFFFFFFu;
    std::uint32_t got_short2    = 0xFFFFFFFFu;
    pfd::popcount_kinds(km.data(), n,
                        &got_partition, &got_short1, &got_short2);
    if (got_partition != expected_partition
        || got_short1 != expected_short1
        || got_short2 != expected_short2) {
        std::fprintf(stderr,
            "FAIL [%s]: n=%u, partition got=%u expected=%u, "
            "short1 got=%u expected=%u, short2 got=%u expected=%u\n",
            label, n,
            got_partition, expected_partition,
            got_short1, expected_short1,
            got_short2, expected_short2);
    }
    CHECK_EQ(got_partition, expected_partition);
    CHECK_EQ(got_short1,    expected_short1);
    CHECK_EQ(got_short2,    expected_short2);
}

void test_empty() {
    std::uint32_t p = 99, s1 = 99, s2 = 99;
    pfd::popcount_kinds(nullptr, 0, &p, &s1, &s2);
    CHECK_EQ(p, 0u);
    CHECK_EQ(s1, 0u);
    CHECK_EQ(s2, 0u);
}

void test_single_each_kind() {
    check_popcount({0}, "single short1");
    check_popcount({1}, "single short2");
    check_popcount({2}, "single partition");
}

void test_homogeneous() {
    for (std::uint32_t n : {std::uint32_t(7), std::uint32_t(31), std::uint32_t(32),
                            std::uint32_t(33), std::uint32_t(63), std::uint32_t(64),
                            std::uint32_t(65), std::uint32_t(127), std::uint32_t(128),
                            std::uint32_t(129), std::uint32_t(255), std::uint32_t(256),
                            std::uint32_t(1024)}) {
        for (std::uint8_t kind : {0, 1, 2}) {
            std::vector<std::uint8_t> kinds(n, kind);
            char buf[64];
            std::snprintf(buf, sizeof(buf), "homogeneous kind=%u n=%u", kind, n);
            check_popcount(kinds, buf);
        }
    }
}

void test_alternating() {
    for (std::uint32_t n : {std::uint32_t(15), std::uint32_t(64), std::uint32_t(96),
                            std::uint32_t(200)}) {
        std::vector<std::uint8_t> kinds(n);
        for (std::uint32_t i = 0; i < n; ++i) {
            kinds[i] = static_cast<std::uint8_t>(i % 3);  // 0,1,2,0,1,2,...
        }
        char buf[64];
        std::snprintf(buf, sizeof(buf), "alternating 0/1/2 n=%u", n);
        check_popcount(kinds, buf);
    }
}

void test_random() {
    std::mt19937 rng(0xCAFE);
    for (std::uint32_t n : {std::uint32_t(33), std::uint32_t(150),
                            std::uint32_t(513), std::uint32_t(4097)}) {
        std::vector<std::uint8_t> kinds(n);
        for (std::uint32_t i = 0; i < n; ++i) {
            kinds[i] = static_cast<std::uint8_t>(rng() % 3);
        }
        char buf[64];
        std::snprintf(buf, sizeof(buf), "random tri-state n=%u", n);
        check_popcount(kinds, buf);
    }
}

// Boundary stress: distinct_count just below / at / above the bulk
// 32-pair chunk boundary to catch off-by-ones in the tail loop.
void test_chunk_boundaries() {
    for (std::uint32_t n : {std::uint32_t(31), std::uint32_t(32), std::uint32_t(33),
                            std::uint32_t(63), std::uint32_t(64), std::uint32_t(65),
                            std::uint32_t(95), std::uint32_t(96), std::uint32_t(97)}) {
        std::vector<std::uint8_t> kinds(n);
        // Mostly partition, sprinkled with one of each other kind to
        // exercise the bulk + tail interaction.
        for (std::uint32_t i = 0; i < n; ++i) kinds[i] = 2;
        if (n >= 1) kinds[0] = 0;
        if (n >= 2) kinds[n - 1] = 1;
        char buf[64];
        std::snprintf(buf, sizeof(buf), "chunk boundary n=%u", n);
        check_popcount(kinds, buf);
    }
}

} // anonymous namespace

int main() {
    test_empty();
    test_single_each_kind();
    test_homogeneous();
    test_alternating();
    test_random();
    test_chunk_boundaries();
    TEST_SUMMARY();
    return g_fail_count == 0 ? 0 : 1;
}
