#include "test_util.hpp"
#include "index/seq_id_dedup.hpp"

#include <cstdint>
#include <cstdio>
#include <vector>

using namespace ikafssn;

static void test_empty_input() {
    std::vector<uint32_t> in;
    std::vector<uint32_t> distinct(0);
    std::vector<uint32_t> occ(0);
    uint32_t d = seq_id_dedup::dedup_seq_ids(in.data(), 0,
                                              distinct.data(), occ.data());
    CHECK_EQ(d, 0u);
}

static void test_single_element() {
    std::vector<uint32_t> in = {42};
    std::vector<uint32_t> distinct(1, 0);
    std::vector<uint32_t> occ(1, 0);
    uint32_t d = seq_id_dedup::dedup_seq_ids(in.data(), 1,
                                              distinct.data(), occ.data());
    CHECK_EQ(d, 1u);
    CHECK_EQ(distinct[0], 42u);
    CHECK_EQ(occ[0], 1u);
}

static void test_all_distinct() {
    std::vector<uint32_t> in = {1, 2, 3, 4, 5};
    std::vector<uint32_t> distinct(5, 0);
    std::vector<uint32_t> occ(5, 0);
    uint32_t d = seq_id_dedup::dedup_seq_ids(in.data(), 5,
                                              distinct.data(), occ.data());
    CHECK_EQ(d, 5u);
    for (uint32_t i = 0; i < 5; i++) {
        CHECK_EQ(distinct[i], i + 1);
        CHECK_EQ(occ[i], 1u);
    }
}

static void test_all_same() {
    std::vector<uint32_t> in(10, 7);
    std::vector<uint32_t> distinct(10, 0);
    std::vector<uint32_t> occ(10, 0);
    uint32_t d = seq_id_dedup::dedup_seq_ids(in.data(), 10,
                                              distinct.data(), occ.data());
    CHECK_EQ(d, 1u);
    CHECK_EQ(distinct[0], 7u);
    CHECK_EQ(occ[0], 10u);
}

static void test_mixed_runs() {
    // {1,1,1,2,3,3,4,4,4,4,5} → distinct {1,2,3,4,5}, occ {3,1,2,4,1}
    std::vector<uint32_t> in = {1, 1, 1, 2, 3, 3, 4, 4, 4, 4, 5};
    std::vector<uint32_t> distinct(in.size(), 0);
    std::vector<uint32_t> occ(in.size(), 0);
    uint32_t d = seq_id_dedup::dedup_seq_ids(in.data(),
                                              static_cast<uint32_t>(in.size()),
                                              distinct.data(), occ.data());
    CHECK_EQ(d, 5u);
    CHECK_EQ(distinct[0], 1u);  CHECK_EQ(occ[0], 3u);
    CHECK_EQ(distinct[1], 2u);  CHECK_EQ(occ[1], 1u);
    CHECK_EQ(distinct[2], 3u);  CHECK_EQ(occ[2], 2u);
    CHECK_EQ(distinct[3], 4u);  CHECK_EQ(occ[3], 4u);
    CHECK_EQ(distinct[4], 5u);  CHECK_EQ(occ[4], 1u);
}

// Boundary lengths: 128, 255 — exercise SIMD-vectorised edge cases.
static void test_boundary_lengths() {
    for (uint32_t n : {128u, 255u, 256u, 511u, 512u, 1023u}) {
        std::vector<uint32_t> in(n);
        for (uint32_t i = 0; i < n; i++) in[i] = i / 2; // run length 2
        std::vector<uint32_t> distinct(n, 0);
        std::vector<uint32_t> occ(n, 0);
        uint32_t d = seq_id_dedup::dedup_seq_ids(in.data(), n,
                                                  distinct.data(), occ.data());
        const uint32_t expected_d = (n + 1) / 2;
        CHECK_EQ(d, expected_d);
        for (uint32_t i = 0; i < d; i++) {
            CHECK_EQ(distinct[i], i);
        }
        // First (n/2) entries have run length 2; the trailing tail is 1 if odd.
        for (uint32_t i = 0; i < d - (n % 2); i++) {
            CHECK_EQ(occ[i], 2u);
        }
        if (n % 2 == 1) {
            CHECK_EQ(occ[d - 1], 1u);
        }
    }
}

// Run length exactly 255.
static void test_run_length_255() {
    std::vector<uint32_t> in;
    for (uint32_t i = 0; i < 255; i++) in.push_back(0);
    in.push_back(1);
    std::vector<uint32_t> distinct(in.size(), 0);
    std::vector<uint32_t> occ(in.size(), 0);
    uint32_t d = seq_id_dedup::dedup_seq_ids(in.data(),
                                              static_cast<uint32_t>(in.size()),
                                              distinct.data(), occ.data());
    CHECK_EQ(d, 2u);
    CHECK_EQ(distinct[0], 0u);
    CHECK_EQ(occ[0], 255u);
    CHECK_EQ(distinct[1], 1u);
    CHECK_EQ(occ[1], 1u);
}

// Regression: u32 occ_count must NOT saturate at 255.  Earlier dedup
// versions silently truncated runs > 255 which corrupted the .kpx
// pos_cursor walk in encode_posting_kpx.
static void test_run_length_above_255() {
    for (uint32_t run : {256u, 1000u, 100000u}) {
        std::vector<uint32_t> in(run, 7u);
        in.push_back(8u);
        std::vector<uint32_t> distinct(in.size(), 0);
        std::vector<uint32_t> occ(in.size(), 0);
        uint32_t d = seq_id_dedup::dedup_seq_ids(in.data(),
                                                  static_cast<uint32_t>(in.size()),
                                                  distinct.data(), occ.data());
        CHECK_EQ(d, 2u);
        CHECK_EQ(distinct[0], 7u);
        CHECK_EQ(occ[0], run);     // must be exact, not saturated
        CHECK_EQ(distinct[1], 8u);
        CHECK_EQ(occ[1], 1u);
    }
}

int main() {
    std::fprintf(stderr, "active dedup tier: %s\n",
                 seq_id_dedup::active_tier_name());
    test_empty_input();
    test_single_element();
    test_all_distinct();
    test_all_same();
    test_mixed_runs();
    test_boundary_lengths();
    test_run_length_255();
    test_run_length_above_255();
    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
