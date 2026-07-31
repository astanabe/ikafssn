#include "test_util.hpp"
#include "search/stage2_chaining.hpp"

#include <algorithm>
#include <random>

using namespace ikafssn;

// chain_hits() writes into a caller-owned vector; the tests below mostly want
// one result list per call, so wrap it.
static std::vector<ChainResult> run_chain(const std::vector<Hit>& hits,
                                          SeqId seq_id,
                                          int span,
                                          bool is_reverse,
                                          const Stage2Config& config) {
    std::vector<ChainResult> out;
    chain_hits(hits, seq_id, span, is_reverse, config, out);
    return out;
}

static void test_single_hit() {
    std::fprintf(stderr, "-- test_single_hit\n");

    std::vector<Hit> hits = {{10, 100}};
    Stage2Config config;
    config.min_nhit_diag = 1; // don't filter
    config.min_score = 1;
    config.max_gap = 100;

    auto result = run_chain(hits, 0, 7, false, config);
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
    config.min_nhit_diag = 1;
    config.min_score = 1;
    config.max_gap = 100;

    auto result = run_chain(hits, 0, 7, false, config);
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
    config.min_nhit_diag = 1;
    config.min_score = 1;
    config.max_gap = 100; // gap_q=40, gap_s=90, diag_diff=50, within 100

    auto result = run_chain(hits, 0, 7, false, config);
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
    config.min_nhit_diag = 1;
    config.min_score = 1;
    config.max_gap = 100;

    config.max_nhit_per_subject_mode = 1;  // strict take-N (no ties)
    auto result = run_chain(hits, 0, 7, false, config);
    CHECK_EQ(result.size(), 1u);
    // Best chain is just one hit (score=1)
    CHECK_EQ(result[0].chainscore, 1u);
}

static void test_min_score_filter() {
    std::fprintf(stderr, "-- test_min_score_filter\n");

    std::vector<Hit> hits = {{10, 100}};
    Stage2Config config;
    config.min_nhit_diag = 1;
    config.min_score = 3; // require at least 3
    config.max_gap = 100;

    auto result = run_chain(hits, 0, 7, false, config);
    CHECK_EQ(result.empty(), true); // filtered out
}

static void test_reverse_strand_flag() {
    std::fprintf(stderr, "-- test_reverse_strand_flag\n");

    std::vector<Hit> hits = {{0, 50}, {7, 57}};
    Stage2Config config;
    config.min_nhit_diag = 1;
    config.min_score = 1;
    config.max_gap = 100;

    auto result = run_chain(hits, 5, 7, true, config);
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
    config.min_nhit_diag = 1;
    config.min_score = 1;
    config.max_gap = 100;

    config.max_nhit_per_subject_mode = 1;  // strict take-N (no ties)
    auto result = run_chain(hits, 0, 7, false, config);
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
    config.min_nhit_diag = 1;
    config.min_score = 1;
    config.max_gap = 100;

    config.max_nhit_per_subject_mode = 1;  // strict take-N (no ties)
    auto result = run_chain(hits, 0, 7, false, config);
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
    config.min_nhit_diag = 1;
    config.min_score = 1;
    config.max_gap = 100;

    config.max_nhit_per_subject_mode = 1;  // strict take-N (no ties)
    auto result = run_chain(hits, 0, 7, false, config);
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
    config.min_nhit_diag = 1;
    config.min_score = 1;
    config.max_gap = 100;
    config.chain_max_lookback = 4;

    auto result = run_chain(hits, 0, 7, false, config);
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
    config.min_nhit_diag = 1;
    config.min_score = 1;
    config.max_gap = 100;

    // B=1: each hit can only look back 1 position.
    // 0: dp=1.  1: pred [0] s_pos 90 > 57, no increase -> dp=1.
    // 2: pred [1] s_pos 57 < 104, gaps ok -> dp=2.  3: pred [2] 104 > 71 -> dp=1.
    // 4: pred [3] 71 < 118, gaps ok -> dp=2.  Best = 2.
    config.max_nhit_per_subject_mode = 1;  // strict take-N (no ties)
    config.chain_max_lookback = 1;
    auto result1 = run_chain(hits, 0, 7, false, config);
    CHECK_EQ(result1.size(), 1u);
    CHECK_EQ(result1[0].chainscore, 2u);

    // B=2: each hit looks back 2 positions.
    // Index 2: looks at [0,1] -> from 0: s 104>90 yes, diag_diff=|104-14-(90-0)|=0<=100. dp=2
    // Index 4: looks at [2,3] -> from 2: s 118>104 yes, diag_diff=|118-28-(104-14)|=0<=100. dp=3
    // Best = 3 (chain A fully recovered)
    config.chain_max_lookback = 2;
    auto result2 = run_chain(hits, 0, 7, false, config);
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
    config.min_nhit_diag = 1;
    config.min_score = 1;
    config.max_gap = 100;
    config.chain_max_lookback = 0; // unlimited

    auto result = run_chain(hits, 0, 7, false, config);
    CHECK_EQ(result.size(), 1u);
    CHECK_EQ(result[0].chainscore, 5u);
}

static void test_empty_hits() {
    std::fprintf(stderr, "-- test_empty_hits\n");

    std::vector<Hit> hits;
    Stage2Config config;
    config.min_nhit_diag = 1;
    config.min_score = 1;
    config.max_gap = 100;

    auto result = run_chain(hits, 0, 7, false, config);
    CHECK_EQ(result.empty(), true);
}

static void test_duplicate_hits_dedup() {
    std::fprintf(stderr, "-- test_duplicate_hits_dedup\n");

    // Duplicate (q_pos, s_pos) pairs from degenerate base expansion
    // should be deduplicated before chaining
    std::vector<Hit> hits = {{10, 100}, {10, 100}, {20, 110}, {20, 110}, {20, 110}};
    Stage2Config config;
    config.min_nhit_diag = 1;
    config.min_score = 1;
    config.max_gap = 100;

    auto result = run_chain(hits, 0, 7, false, config);
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
    config.min_nhit_diag = 1;
    config.min_score = 1;
    config.max_gap = 50; // tight gap: prevents cross-region chaining
    config.max_nhit_per_subject = 2;

    auto result = run_chain(hits, 0, 7, false, config);
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
    config.min_nhit_diag = 1;
    config.min_score = 1;
    config.max_gap = 50;
    config.max_nhit_per_subject = 0; // unlimited

    auto result = run_chain(hits, 0, 7, false, config);
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
    config.min_nhit_diag = 1;
    config.min_score = 1;
    config.max_gap = 50;
    // max_nhit_per_subject = 1 (default)

    auto result = run_chain(hits, 0, 7, false, config);
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
    config.min_nhit_diag = 1;
    config.min_score = 1;
    config.max_gap = 50;
    config.max_nhit_per_subject = 3;

    auto result = run_chain(hits, 0, 7, false, config);
    CHECK_EQ(result.size(), 3u);
    // Greedy removal: first finds best (4), then best of remainder (3), then (2)
    CHECK_EQ(result[0].chainscore, 4u);
    CHECK_EQ(result[1].chainscore, 3u);
    CHECK_EQ(result[2].chainscore, 2u);
}

static void test_tie_inclusive_take_n() {
    std::fprintf(stderr, "-- test_tie_inclusive_take_n\n");

    // Three independent regions: A and B tie at score 3, C is score 2.
    std::vector<Hit> hits = {
        {0, 100}, {7, 107}, {14, 114},    // region A: score 3
        {0, 500}, {7, 507}, {14, 514},    // region B: score 3
        {0, 900}, {7, 907},               // region C: score 2
    };
    Stage2Config config;
    config.min_nhit_diag = 1;
    config.min_score = 1;
    config.max_gap = 50;
    config.max_nhit_per_subject = 1;

    // Strict take-1 (mode 1): only the single best chain survives.
    config.max_nhit_per_subject_mode = 1;
    auto strict = run_chain(hits, 0, 7, false, config);
    CHECK_EQ(strict.size(), 1u);
    CHECK_EQ(strict[0].chainscore, 3u);

    // Tie-inclusive take-1 (mode 3): the top score (3) plus its tie -> 2.
    config.max_nhit_per_subject_mode = 3;
    auto tie = run_chain(hits, 0, 7, false, config);
    CHECK_EQ(tie.size(), 2u);
    CHECK_EQ(tie[0].chainscore, 3u);
    CHECK_EQ(tie[1].chainscore, 3u);

    // Tie-inclusive N=2: 2nd score is 3, so the score-2 region C is dropped.
    config.max_nhit_per_subject = 2;
    auto tie2 = run_chain(hits, 0, 7, false, config);
    CHECK_EQ(tie2.size(), 2u);

    // Unlimited (N=0) returns everything regardless of mode.
    config.max_nhit_per_subject = 0;
    auto all = run_chain(hits, 0, 7, false, config);
    CHECK_EQ(all.size(), 3u);
}

// Two collinear regions plus a scatter of off-diagonal noise, enough hits to
// exercise the lookback window and several extraction iterations.
static std::vector<Hit> make_two_region_hits() {
    std::vector<Hit> hits;
    for (uint32_t i = 0; i < 40; i++) hits.push_back({i * 7, 1000 + i * 7});
    for (uint32_t i = 0; i < 25; i++) hits.push_back({i * 7, 5000 + i * 8});
    for (uint32_t i = 0; i < 10; i++) hits.push_back({i * 31 + 3, 9000 - i * 211});
    std::sort(hits.begin(), hits.end(), [](const Hit& a, const Hit& b) {
        return a.q_pos < b.q_pos || (a.q_pos == b.q_pos && a.s_pos < b.s_pos);
    });
    return hits;
}

static Stage2Config two_region_config() {
    Stage2Config config;
    config.min_nhit_diag = 1;
    config.min_score = 1;
    config.max_gap = 50;
    config.max_nhit_per_subject = 0;  // unlimited: exercise the removal loop
    return config;
}

static bool same_chains(const std::vector<ChainResult>& a,
                        const std::vector<ChainResult>& b) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); i++) {
        if (a[i].seq_id != b[i].seq_id) return false;
        if (a[i].chainscore != b[i].chainscore) return false;
        if (a[i].q_start != b[i].q_start || a[i].q_end != b[i].q_end) return false;
        if (a[i].s_start != b[i].s_start || a[i].s_end != b[i].s_end) return false;
        if (a[i].is_reverse != b[i].is_reverse) return false;
    }
    return true;
}

static void test_scratch_reuse_across_calls() {
    std::fprintf(stderr, "-- test_scratch_reuse_across_calls\n");

    // A large subject, then a small one, then the large one again.  The
    // per-thread working buffers are sized for the first call, so a missing
    // clear or resize would let its leftovers leak into the later calls.
    const std::vector<Hit> big = make_two_region_hits();
    const std::vector<Hit> small = {{0, 100}, {7, 107}, {14, 114}};
    const Stage2Config config = two_region_config();

    std::vector<ChainResult> out;
    chain_hits(big, 3, 7, false, config, out);
    const std::vector<ChainResult> first = out;
    CHECK_EQ(first.empty(), false);

    chain_hits(small, 3, 7, false, config, out);
    CHECK_EQ(out.size(), 1u);
    CHECK_EQ(out[0].chainscore, 3u);

    chain_hits(big, 3, 7, false, config, out);
    CHECK_EQ(same_chains(first, out), true);
}

static void test_unsorted_input_equivalence() {
    std::fprintf(stderr, "-- test_unsorted_input_equivalence\n");

    // The same hit set in five arrangements: one ascending run, one descending
    // run, two ascending runs, two descending runs, and a shuffle.
    const std::vector<Hit> ascending = make_two_region_hits();
    const size_t n = ascending.size();

    std::vector<Hit> descending(ascending.rbegin(), ascending.rend());

    std::vector<Hit> two_asc_runs, two_desc_runs;
    for (size_t i = 0; i < n; i += 2) two_asc_runs.push_back(ascending[i]);
    const size_t run1_len = two_asc_runs.size();
    for (size_t i = 1; i < n; i += 2) two_asc_runs.push_back(ascending[i]);
    two_desc_runs = two_asc_runs;
    std::reverse(two_desc_runs.begin(), two_desc_runs.begin() + run1_len);
    std::reverse(two_desc_runs.begin() + run1_len, two_desc_runs.end());

    std::vector<Hit> shuffled = ascending;
    std::mt19937 rng(20260727u);
    std::shuffle(shuffled.begin(), shuffled.end(), rng);

    const Stage2Config config = two_region_config();
    const std::vector<ChainResult> expected =
        run_chain(ascending, 3, 7, false, config);
    CHECK_EQ(expected.empty(), false);

    CHECK_EQ(same_chains(expected, run_chain(descending, 3, 7, false, config)), true);
    CHECK_EQ(same_chains(expected, run_chain(two_asc_runs, 3, 7, false, config)), true);
    CHECK_EQ(same_chains(expected, run_chain(two_desc_runs, 3, 7, false, config)), true);
    CHECK_EQ(same_chains(expected, run_chain(shuffled, 3, 7, false, config)), true);
}

static void test_out_is_cleared() {
    std::fprintf(stderr, "-- test_out_is_cleared\n");

    // Stage 2B hands the same vector to every candidate, so whatever the
    // previous subject left behind must not survive into the next call --
    // including on the two early-return paths.
    Stage2Config config;
    config.min_nhit_diag = 1;
    config.min_score = 1;
    config.max_gap = 100;

    ChainResult stale{};
    stale.seq_id = 99;
    stale.chainscore = 12345;

    std::vector<ChainResult> out(3, stale);
    chain_hits({}, 0, 7, false, config, out);      // no hits at all
    CHECK_EQ(out.empty(), true);

    out.assign(3, stale);
    config.min_score = 5;                          // no chain reaches it
    chain_hits({{10, 100}}, 0, 7, false, config, out);
    CHECK_EQ(out.empty(), true);

    out.assign(3, stale);
    config.min_score = 1;
    chain_hits({{10, 100}}, 0, 7, false, config, out);
    CHECK_EQ(out.size(), 1u);
    CHECK_EQ(out[0].seq_id, 0u);
    CHECK_EQ(out[0].chainscore, 1u);
}

static void test_max_gap_boundary() {
    std::fprintf(stderr, "-- test_max_gap_boundary\n");

    // Anchor on diagonal 1000; the partner's diagonal is offset by exactly
    // +/- max_gap and by one more.  Both signs must behave the same way.
    Stage2Config config;
    config.min_nhit_diag = 1;
    config.min_score = 1;
    config.max_gap = 100;
    config.max_nhit_per_subject_mode = 1;  // strict take-N (no ties)

    auto score_for = [&](uint32_t s2) {
        std::vector<Hit> hits = {{0, 1000}, {200, s2}};
        auto result = run_chain(hits, 0, 7, false, config);
        CHECK_EQ(result.size(), 1u);
        return result.empty() ? 0u : result[0].chainscore;
    };

    CHECK_EQ(score_for(1300), 2u);  // diag_diff = +100 == max_gap
    CHECK_EQ(score_for(1301), 1u);  // diag_diff = +101 >  max_gap
    CHECK_EQ(score_for(1100), 2u);  // diag_diff = -100, |.| == max_gap
    CHECK_EQ(score_for(1099), 1u);  // diag_diff = -101, |.| >  max_gap
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
    test_tie_inclusive_take_n();
    test_scratch_reuse_across_calls();
    test_unsorted_input_equivalence();
    test_out_is_cleared();
    test_max_gap_boundary();

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
