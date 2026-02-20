#include "test_util.hpp"
#include "search/diagonal_filter.hpp"

using namespace ikafssn;

static void test_no_filter_when_threshold_1() {
    std::fprintf(stderr, "-- test_no_filter_when_threshold_1\n");

    std::vector<Hit> hits = {
        {0, 100}, {5, 105}, {10, 200}, {15, 300}
    };

    auto result = diagonal_filter(hits, 1);
    CHECK_EQ(result.size(), hits.size());
}

static void test_filter_isolates() {
    std::fprintf(stderr, "-- test_filter_isolates\n");

    // Hits on diagonal 100 (s-q=100): {0,100}, {5,105}, {10,110}
    // Hits on diagonal 200 (s-q=200): {10,210} -- only 1 hit
    std::vector<Hit> hits = {
        {0, 100}, {5, 105}, {10, 110}, {10, 210}
    };

    auto result = diagonal_filter(hits, 2);
    CHECK_EQ(result.size(), 3u); // only diagonal 100 hits survive
    for (const auto& h : result) {
        int32_t diag = static_cast<int32_t>(h.s_pos) - static_cast<int32_t>(h.q_pos);
        CHECK_EQ(diag, 100);
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

    auto result = diagonal_filter(hits, 3);
    CHECK_EQ(result.size(), 3u); // only diagonal 50 hits survive
}

static void test_empty_input() {
    std::fprintf(stderr, "-- test_empty_input\n");

    std::vector<Hit> hits;
    auto result = diagonal_filter(hits, 2);
    CHECK_EQ(result.size(), 0u);
}

static void test_negative_diagonal() {
    std::fprintf(stderr, "-- test_negative_diagonal\n");

    // Diagonal -10 (s_pos < q_pos): {20, 10}, {25, 15}
    std::vector<Hit> hits = {
        {20, 10}, {25, 15}, {30, 100}
    };

    auto result = diagonal_filter(hits, 2);
    CHECK_EQ(result.size(), 2u);
    for (const auto& h : result) {
        int32_t diag = static_cast<int32_t>(h.s_pos) - static_cast<int32_t>(h.q_pos);
        CHECK_EQ(diag, -10);
    }
}

int main() {
    test_no_filter_when_threshold_1();
    test_filter_isolates();
    test_filter_with_higher_threshold();
    test_empty_input();
    test_negative_diagonal();

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
