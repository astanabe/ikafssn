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

    auto result = chain_hits(hits, 0, 7, false, config);
    CHECK_EQ(result.size(), 1u);
    CHECK_EQ(result[0].chainscore, 1u);
    CHECK_EQ(result[0].q_start, 10u);
    CHECK_EQ(result[0].q_end, 17u); // 10 + k
    CHECK_EQ(result[0].s_start, 100u);
    CHECK_EQ(result[0].s_end, 107u);
    CHECK_EQ(result[0].is_reverse, false);
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

    auto result = chain_hits(hits, 0, 7, false, config);
    CHECK_EQ(result.size(), 1u);
    CHECK_EQ(result[0].chainscore, 5u);
    CHECK_EQ(result[0].q_start, 0u);
    CHECK_EQ(result[0].q_end, 35u); // 28 + 7
    CHECK_EQ(result[0].s_start, 90u);
    CHECK_EQ(result[0].s_end, 125u); // 118 + 7
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

    auto result = chain_hits(hits, 0, 7, false, config);
    CHECK_EQ(result.size(), 1u);
    CHECK_EQ(result[0].chainscore, 3u);
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

    auto result = chain_hits(hits, 0, 7, false, config);
    CHECK_EQ(result.size(), 1u);
    // Best chain is just one hit (score=1)
    CHECK_EQ(result[0].chainscore, 1u);
}

static void test_min_score_filter() {
    std::fprintf(stderr, "-- test_min_score_filter\n");

    std::vector<Hit> hits = {{10, 100}};
    Stage2Config config;
    config.min_diag_hits = 1;
    config.min_score = 3; // require at least 3
    config.max_gap = 100;

    auto result = chain_hits(hits, 0, 7, false, config);
    CHECK_EQ(result.empty(), true); // filtered out
}

static void test_reverse_strand_flag() {
    std::fprintf(stderr, "-- test_reverse_strand_flag\n");

    std::vector<Hit> hits = {{0, 50}, {7, 57}};
    Stage2Config config;
    config.min_diag_hits = 1;
    config.min_score = 1;
    config.max_gap = 100;

    auto result = chain_hits(hits, 5, 7, true, config);
    CHECK_EQ(result.size(), 1u);
    CHECK_EQ(result[0].is_reverse, true);
    CHECK_EQ(result[0].seq_id, 5u);
    CHECK_EQ(result[0].chainscore, 2u);
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

    auto result = chain_hits(hits, 0, 7, false, config);
    CHECK_EQ(result.size(), 1u);
    // Each hit is independent, best chain score = 1
    CHECK_EQ(result[0].chainscore, 1u);
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

    auto result = chain_hits(hits, 0, 7, false, config);
    CHECK_EQ(result.size(), 1u);
    // Only one distinct q_pos, so the longest chain must be 1
    CHECK_EQ(result[0].chainscore, 1u);
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

    auto result = chain_hits(hits, 0, 7, false, config);
    CHECK_EQ(result.size(), 1u);
    // Best chain: pick one hit from q_pos=0 and one from q_pos=20 → score 2
    CHECK_EQ(result[0].chainscore, 2u);
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

    auto result = chain_hits(hits, 0, 7, false, config);
    CHECK_EQ(result.size(), 1u);
    CHECK_EQ(result[0].chainscore, 5u);
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
    auto result1 = chain_hits(hits, 0, 7, false, config);
    CHECK_EQ(result1.size(), 1u);
    CHECK_EQ(result1[0].chainscore, 2u);

    // B=2: each hit looks back 2 positions.
    // Index 2: looks at [0,1] -> from 0: s 104>90 yes, diag_diff=|104-14-(90-0)|=0<=100. dp=2
    // Index 4: looks at [2,3] -> from 2: s 118>104 yes, diag_diff=|118-28-(104-14)|=0<=100. dp=3
    // Best = 3 (chain A fully recovered)
    config.chain_max_lookback = 2;
    auto result2 = chain_hits(hits, 0, 7, false, config);
    CHECK_EQ(result2.size(), 1u);
    CHECK_EQ(result2[0].chainscore, 3u);
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

    auto result = chain_hits(hits, 0, 7, false, config);
    CHECK_EQ(result.size(), 1u);
    CHECK_EQ(result[0].chainscore, 5u);
}

static void test_empty_hits() {
    std::fprintf(stderr, "-- test_empty_hits\n");

    std::vector<Hit> hits;
    Stage2Config config;
    config.min_diag_hits = 1;
    config.min_score = 1;
    config.max_gap = 100;

    auto result = chain_hits(hits, 0, 7, false, config);
    CHECK_EQ(result.empty(), true);
}

static void test_duplicate_hits_dedup() {
    std::fprintf(stderr, "-- test_duplicate_hits_dedup\n");

    // Duplicate (q_pos, s_pos) pairs from degenerate base expansion
    // should be deduplicated before chaining
    std::vector<Hit> hits = {{10, 100}, {10, 100}, {20, 110}, {20, 110}, {20, 110}};
    Stage2Config config;
    config.min_diag_hits = 1;
    config.min_score = 1;
    config.max_gap = 100;

    auto result = chain_hits(hits, 0, 7, false, config);
    CHECK_EQ(result.size(), 1u);
    CHECK_EQ(result[0].chainscore, 2u);  // 2 distinct positions, not 5
}

static void test_multi_chain_two_regions() {
    std::fprintf(stderr, "-- test_multi_chain_two_regions\n");

    // Two independent collinear regions on the same subject
    // Region A: q_pos 0-28, s_pos 90-118 (diag 90)
    // Region B: q_pos 0-28, s_pos 500-528 (diag 500)
    // These share q_pos values but have very different s_pos, so gap constraint
    // prevents them from being in the same chain.
    std::vector<Hit> hits = {
        {0, 90}, {7, 97}, {14, 104}, {21, 111}, {28, 118},
        {0, 500}, {7, 507}, {14, 514}, {21, 521}, {28, 528}
    };
    Stage2Config config;
    config.min_diag_hits = 1;
    config.min_score = 1;
    config.max_gap = 50; // tight gap: prevents cross-region chaining
    config.max_nhit_per_subject = 2;

    auto result = chain_hits(hits, 0, 7, false, config);
    CHECK_EQ(result.size(), 2u);
    CHECK_EQ(result[0].chainscore, 5u);
    CHECK_EQ(result[1].chainscore, 5u);
}

static void test_multi_chain_unlimited() {
    std::fprintf(stderr, "-- test_multi_chain_unlimited\n");

    // Three independent regions
    std::vector<Hit> hits = {
        {0, 100}, {7, 107}, {14, 114},    // region A: score 3
        {0, 500}, {7, 507},                // region B: score 2
        {0, 900}, {7, 907}, {14, 914}, {21, 921} // region C: score 4
    };
    Stage2Config config;
    config.min_diag_hits = 1;
    config.min_score = 1;
    config.max_gap = 50;
    config.max_nhit_per_subject = 0; // unlimited

    auto result = chain_hits(hits, 0, 7, false, config);
    CHECK_EQ(result.size(), 3u);
    // Should be in descending score order: 4, 3, 2
    CHECK_EQ(result[0].chainscore, 4u);
    CHECK_EQ(result[1].chainscore, 3u);
    CHECK_EQ(result[2].chainscore, 2u);
}

static void test_multi_chain_default_one() {
    std::fprintf(stderr, "-- test_multi_chain_default_one\n");

    // With default max_nhit_per_subject=1, only the best chain is returned
    std::vector<Hit> hits = {
        {0, 100}, {7, 107}, {14, 114},
        {0, 500}, {7, 507}, {14, 514}, {21, 521}
    };
    Stage2Config config;
    config.min_diag_hits = 1;
    config.min_score = 1;
    config.max_gap = 50;
    // max_nhit_per_subject = 1 (default)

    auto result = chain_hits(hits, 0, 7, false, config);
    CHECK_EQ(result.size(), 1u);
    CHECK_EQ(result[0].chainscore, 4u); // best chain
}

static void test_multi_chain_score_order() {
    std::fprintf(stderr, "-- test_multi_chain_score_order\n");

    // Verify chains are returned in score-descending order (greedy best-first)
    std::vector<Hit> hits = {
        {0, 100}, {7, 107},                       // region A: score 2
        {0, 500}, {7, 507}, {14, 514}, {21, 521}, // region B: score 4
        {0, 900}, {7, 907}, {14, 914}              // region C: score 3
    };
    Stage2Config config;
    config.min_diag_hits = 1;
    config.min_score = 1;
    config.max_gap = 50;
    config.max_nhit_per_subject = 3;

    auto result = chain_hits(hits, 0, 7, false, config);
    CHECK_EQ(result.size(), 3u);
    // Greedy removal: first finds best (4), then best of remainder (3), then (2)
    CHECK_EQ(result[0].chainscore, 4u);
    CHECK_EQ(result[1].chainscore, 3u);
    CHECK_EQ(result[2].chainscore, 2u);
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
    test_duplicate_hits_dedup();
    test_multi_chain_two_regions();
    test_multi_chain_unlimited();
    test_multi_chain_default_one();
    test_multi_chain_score_order();

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
