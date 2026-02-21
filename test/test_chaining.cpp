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

static void test_same_qpos_not_chained() {
    std::fprintf(stderr, "-- test_same_qpos_not_chained\n");

    // Multiple hits at the same q_pos (a single query k-mer matching
    // several subject positions). These must NOT inflate the chain score.
    std::vector<Hit> hits = {
        {10, 100}, {10, 110}, {10, 120}, {10, 130}
    };
    Stage2Config config;
    config.min_diag_hits = 1;
    config.min_score = 1;
    config.max_gap = 100;

    ChainResult cr = chain_hits(hits, 0, 7, false, config);
    // Only one distinct q_pos, so the longest chain must be 1
    CHECK_EQ(cr.score, 1u);
}

static void test_same_qpos_mixed_with_distinct() {
    std::fprintf(stderr, "-- test_same_qpos_mixed_with_distinct\n");

    // Two distinct q_pos values (0 and 20), each with multiple s_pos hits.
    // The chain should use at most one hit per q_pos.
    std::vector<Hit> hits = {
        {0, 100}, {0, 110}, {0, 120},
        {20, 200}, {20, 210}, {20, 220}
    };
    Stage2Config config;
    config.min_diag_hits = 1;
    config.min_score = 1;
    config.max_gap = 100;

    ChainResult cr = chain_hits(hits, 0, 7, false, config);
    // Best chain: pick one hit from q_pos=0 and one from q_pos=20 → score 2
    CHECK_EQ(cr.score, 2u);
}

static void test_chain_max_lookback_basic() {
    std::fprintf(stderr, "-- test_chain_max_lookback_basic\n");

    // 5 collinear hits on diagonal 90, all within lookback window B=4
    std::vector<Hit> hits = {
        {0, 90}, {7, 97}, {14, 104}, {21, 111}, {28, 118}
    };
    Stage2Config config;
    config.min_diag_hits = 1;
    config.min_score = 1;
    config.max_gap = 100;
    config.chain_max_lookback = 4;

    ChainResult cr = chain_hits(hits, 0, 7, false, config);
    CHECK_EQ(cr.score, 5u);
}

static void test_chain_max_lookback_interleaved() {
    std::fprintf(stderr, "-- test_chain_max_lookback_interleaved\n");

    // Two interleaved chains on different diagonals:
    // Chain A (diag 90): {0,90}, {14,104}, {28,118}
    // Chain B (diag 50): {7,57}, {21,71}
    // Sorted by q_pos: {0,90}, {7,57}, {14,104}, {21,71}, {28,118}
    std::vector<Hit> hits = {
        {0, 90}, {7, 57}, {14, 104}, {21, 71}, {28, 118}
    };
    Stage2Config config;
    config.min_diag_hits = 1;
    config.min_score = 1;
    config.max_gap = 100;

    // B=1: each hit can only look back 1 position.
    // Index 0: dp=1
    // Index 1: looks at [0] only -> diag diff |57-7-(90-0)|=|50-90|=40 <=100, but s needs increase: 57<90? No 57<90 -> s_pos[1]=57 < s_pos[0]=90, skip. dp=1
    // Index 2: looks at [1] only -> s_pos[2]=104>57 yes, gap check: gap_q=14-7=7, gap_s=104-57=47, diag_diff=40<=100. dp=2
    // Index 3: looks at [2] only -> s_pos[3]=71 < s_pos[2]=104, skip. dp=1
    // Index 4: looks at [3] only -> s_pos[4]=118>71 yes, gap_q=28-21=7, gap_s=118-71=47, diag_diff=40<=100. dp=2
    // Best = 2
    config.chain_max_lookback = 1;
    ChainResult cr1 = chain_hits(hits, 0, 7, false, config);
    CHECK_EQ(cr1.score, 2u);

    // B=2: each hit looks back 2 positions.
    // Index 2: looks at [0,1] -> from 0: s 104>90 yes, diag_diff=|104-14-(90-0)|=0<=100. dp=2
    // Index 4: looks at [2,3] -> from 2: s 118>104 yes, diag_diff=|118-28-(104-14)|=0<=100. dp=3
    // Best = 3 (chain A fully recovered)
    config.chain_max_lookback = 2;
    ChainResult cr2 = chain_hits(hits, 0, 7, false, config);
    CHECK_EQ(cr2.score, 3u);
}

static void test_chain_max_lookback_zero_unlimited() {
    std::fprintf(stderr, "-- test_chain_max_lookback_zero_unlimited\n");

    // Same as test_perfect_chain but with chain_max_lookback=0 (unlimited)
    std::vector<Hit> hits = {
        {0, 90}, {7, 97}, {14, 104}, {21, 111}, {28, 118}
    };
    Stage2Config config;
    config.min_diag_hits = 1;
    config.min_score = 1;
    config.max_gap = 100;
    config.chain_max_lookback = 0; // unlimited

    ChainResult cr = chain_hits(hits, 0, 7, false, config);
    CHECK_EQ(cr.score, 5u);
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
    test_same_qpos_not_chained();
    test_same_qpos_mixed_with_distinct();
    test_chain_max_lookback_basic();
    test_chain_max_lookback_interleaved();
    test_chain_max_lookback_zero_unlimited();
    test_empty_hits();

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
