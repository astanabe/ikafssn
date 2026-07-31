#include "test_util.hpp"
#include "search/diagonal_filter.hpp"

using namespace ikafssn;

// diagonal_filter() shrinks its input in place and borrows a working buffer;
// the tests below mostly want one filtered list per call, so wrap it.
static std::vector<Hit> run_filter(std::vector<Hit> hits, uint32_t min_nhit_diag) {
    std::vector<int32_t> diag_scratch;
    diagonal_filter(hits, min_nhit_diag, diag_scratch);
    return hits;
}

static int32_t diag_of(const Hit& h) {
    return static_cast<int32_t>(h.s_pos) - static_cast<int32_t>(h.q_pos);
}

static void test_no_filter_when_threshold_1() {
    std::fprintf(stderr, "-- test_no_filter_when_threshold_1\n");

    std::vector<Hit> hits = {
        {0, 100}, {5, 105}, {10, 200}, {15, 300}
    };

    auto result = run_filter(hits, 1);
    CHECK_EQ(result.size(), hits.size());
}

static void test_filter_isolates() {
    std::fprintf(stderr, "-- test_filter_isolates\n");

    // Hits on diagonal 100 (s-q=100): {0,100}, {5,105}, {10,110}
    // Hits on diagonal 200 (s-q=200): {10,210} -- only 1 hit
    std::vector<Hit> hits = {
        {0, 100}, {5, 105}, {10, 110}, {10, 210}
    };

    auto result = run_filter(hits, 2);
    CHECK_EQ(result.size(), 3u); // only diagonal 100 hits survive
    for (const auto& h : result) {
        CHECK_EQ(diag_of(h), 100);
    }
}

static void test_filter_with_higher_threshold() {
    std::fprintf(stderr, "-- test_filter_with_higher_threshold\n");

    // Diagonal 50: 3 hits
    // Diagonal 100: 2 hits
    std::vector<Hit> hits = {
        {0, 50}, {5, 55}, {10, 60},  // diag=50
        {0, 100}, {5, 105},           // diag=100
    };

    auto result = run_filter(hits, 3);
    CHECK_EQ(result.size(), 3u); // only diagonal 50 hits survive
}

static void test_empty_input() {
    std::fprintf(stderr, "-- test_empty_input\n");

    std::vector<Hit> hits;
    auto result = run_filter(hits, 2);
    CHECK_EQ(result.size(), 0u);
}

static void test_negative_diagonal() {
    std::fprintf(stderr, "-- test_negative_diagonal\n");

    // Diagonal -10 (s_pos < q_pos): {20, 10}, {25, 15}
    std::vector<Hit> hits = {
        {20, 10}, {25, 15}, {30, 100}
    };

    auto result = run_filter(hits, 2);
    CHECK_EQ(result.size(), 2u);
    for (const auto& h : result) {
        CHECK_EQ(diag_of(h), -10);
    }
}

static void test_input_order_preserved() {
    std::fprintf(stderr, "-- test_input_order_preserved\n");

    // Survivors are interleaved across three diagonals so any reordering by
    // diagonal would show up.  chain_hits() feeds a (q_pos, s_pos)-ordered
    // list and relies on the filter leaving it ordered.
    std::vector<Hit> hits = {
        {0, 100},   // diag 100, kept
        {0, 500},   // diag 500, dropped (alone)
        {1, 51},    // diag 50,  kept
        {2, 102},   // diag 100, kept
        {3, 53},    // diag 50,  kept
        {4, 900},   // diag 896, dropped (alone)
        {5, 105},   // diag 100, kept
    };

    auto result = run_filter(hits, 2);
    CHECK_EQ(result.size(), 5u);
    const std::vector<Hit> expected = {
        {0, 100}, {1, 51}, {2, 102}, {3, 53}, {5, 105}
    };
    for (size_t i = 0; i < expected.size() && i < result.size(); i++) {
        CHECK_EQ(result[i].q_pos, expected[i].q_pos);
        CHECK_EQ(result[i].s_pos, expected[i].s_pos);
    }
}

static void test_scratch_reuse_across_calls() {
    std::fprintf(stderr, "-- test_scratch_reuse_across_calls\n");

    // A long call followed by a short one must not let the first call's
    // leftover diagonals inflate the second call's counts.
    std::vector<Hit> big;
    for (uint32_t i = 0; i < 64; i++) big.push_back({i, i + 7});  // diag 7
    std::vector<Hit> small = { {0, 7}, {10, 200} };               // diag 7, 190

    std::vector<int32_t> diag_scratch;

    std::vector<Hit> a = big;
    diagonal_filter(a, 2, diag_scratch);
    CHECK_EQ(a.size(), 64u);

    std::vector<Hit> b = small;
    diagonal_filter(b, 2, diag_scratch);
    CHECK_EQ(b.size(), 0u);  // neither diagonal reaches 2 hits on its own

    std::vector<Hit> c = big;
    diagonal_filter(c, 2, diag_scratch);
    CHECK_EQ(c.size(), 64u);

    // Same input, fresh buffer -> same output.
    auto fresh = run_filter(big, 2);
    CHECK_EQ(c.size(), fresh.size());
    for (size_t i = 0; i < fresh.size() && i < c.size(); i++) {
        CHECK_EQ(c[i].q_pos, fresh[i].q_pos);
        CHECK_EQ(c[i].s_pos, fresh[i].s_pos);
    }
}

int main() {
    test_no_filter_when_threshold_1();
    test_filter_isolates();
    test_filter_with_higher_threshold();
    test_empty_input();
    test_negative_diagonal();
    test_input_order_preserved();
    test_scratch_reuse_across_calls();

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
