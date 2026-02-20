#include "test_util.hpp"
#include "search/stage2_chaining.hpp"

using namespace ikafssn;

static void test_single_hit() {
    std::fprintf(stderr, "-- test_single_hit\n");

    std::vector<Hit> hits = {{10, 100}};
    Stage2Config config;
    config.min_diag_hits = 1; // don't filter
    config.min_score = 1;
    config.max_gap = 100;

    ChainResult cr = chain_hits(hits, 0, 7, false, config);
    CHECK_EQ(cr.score, 1u);
    CHECK_EQ(cr.q_start, 10u);
    CHECK_EQ(cr.q_end, 17u); // 10 + k
    CHECK_EQ(cr.s_start, 100u);
    CHECK_EQ(cr.s_end, 107u);
    CHECK_EQ(cr.is_reverse, false);
}

static void test_perfect_chain() {
    std::fprintf(stderr, "-- test_perfect_chain\n");

    // 5 collinear hits on diagonal 90 (s-q=90), spaced k=7 apart
    std::vector<Hit> hits = {
        {0, 90}, {7, 97}, {14, 104}, {21, 111}, {28, 118}
    };
    Stage2Config config;
    config.min_diag_hits = 1;
    config.min_score = 1;
    config.max_gap = 100;

    ChainResult cr = chain_hits(hits, 0, 7, false, config);
    CHECK_EQ(cr.score, 5u);
    CHECK_EQ(cr.q_start, 0u);
    CHECK_EQ(cr.q_end, 35u); // 28 + 7
    CHECK_EQ(cr.s_start, 90u);
    CHECK_EQ(cr.s_end, 125u); // 118 + 7
}

static void test_chain_with_gap() {
    std::fprintf(stderr, "-- test_chain_with_gap\n");

    // 3 hits: first two on same diagonal, third shifted but within max_gap
    std::vector<Hit> hits = {
        {0, 100}, {10, 110}, {50, 200}
    };
    Stage2Config config;
    config.min_diag_hits = 1;
    config.min_score = 1;
    config.max_gap = 100; // gap_q=40, gap_s=90, diag_diff=50, within 100

    ChainResult cr = chain_hits(hits, 0, 7, false, config);
    CHECK_EQ(cr.score, 3u);
}

static void test_chain_gap_exceeded() {
    std::fprintf(stderr, "-- test_chain_gap_exceeded\n");

    // Two hits: diagonal difference exceeds max_gap
    std::vector<Hit> hits = {
        {0, 100}, {10, 300} // gap_q=10, gap_s=200, diag_diff=190 > max_gap=100
    };
    Stage2Config config;
    config.min_diag_hits = 1;
    config.min_score = 1;
    config.max_gap = 100;

    ChainResult cr = chain_hits(hits, 0, 7, false, config);
    // Best chain is just one hit (score=1)
    CHECK_EQ(cr.score, 1u);
}

static void test_min_score_filter() {
    std::fprintf(stderr, "-- test_min_score_filter\n");

    std::vector<Hit> hits = {{10, 100}};
    Stage2Config config;
    config.min_diag_hits = 1;
    config.min_score = 3; // require at least 3
    config.max_gap = 100;

    ChainResult cr = chain_hits(hits, 0, 7, false, config);
    CHECK_EQ(cr.score, 0u); // filtered out
}

static void test_reverse_strand_flag() {
    std::fprintf(stderr, "-- test_reverse_strand_flag\n");

    std::vector<Hit> hits = {{0, 50}, {7, 57}};
    Stage2Config config;
    config.min_diag_hits = 1;
    config.min_score = 1;
    config.max_gap = 100;

    ChainResult cr = chain_hits(hits, 5, 7, true, config);
    CHECK_EQ(cr.is_reverse, true);
    CHECK_EQ(cr.seq_id, 5u);
    CHECK_EQ(cr.score, 2u);
}

static void test_non_collinear_hits() {
    std::fprintf(stderr, "-- test_non_collinear_hits\n");

    // Hits where s_pos decreases as q_pos increases (not collinear)
    std::vector<Hit> hits = {
        {0, 200}, {10, 100}, {20, 50}
    };
    Stage2Config config;
    config.min_diag_hits = 1;
    config.min_score = 1;
    config.max_gap = 100;

    ChainResult cr = chain_hits(hits, 0, 7, false, config);
    // Each hit is independent, best chain score = 1
    CHECK_EQ(cr.score, 1u);
}

static void test_empty_hits() {
    std::fprintf(stderr, "-- test_empty_hits\n");

    std::vector<Hit> hits;
    Stage2Config config;
    config.min_diag_hits = 1;
    config.min_score = 1;
    config.max_gap = 100;

    ChainResult cr = chain_hits(hits, 0, 7, false, config);
    CHECK_EQ(cr.score, 0u);
}

int main() {
    test_single_hit();
    test_perfect_chain();
    test_chain_with_gap();
    test_chain_gap_exceeded();
    test_min_score_filter();
    test_reverse_strand_flag();
    test_non_collinear_hits();
    test_empty_hits();

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
