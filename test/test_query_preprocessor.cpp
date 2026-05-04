#include "test_util.hpp"
#include "search/query_preprocessor.hpp"
#include "search/volume_searcher.hpp"
#include "protocol/messages.hpp"

#include <cstdio>
#include <string>

using namespace ikafssn;

// Phase 1 spec: Nqkmer = seq_len - span + 1 (pure window count, content-
// independent). When seq_len < span, the query is too short and a single
// k-mer cannot fit in the window — kSkipQueryTooShort fires.
static void test_too_short() {
    std::fprintf(stderr, "-- test_too_short\n");
    SearchConfig config;
    config.stage1.min_stage1_score = 1;

    // k=9, seq_len=8 (< 9): too short
    auto q = preprocess_query<uint16_t>(std::string(8, 'A'), 9, nullptr, config);
    CHECK_EQ(q.skip_reason, (uint8_t)kSkipQueryTooShort);
    CHECK(!q.skip_detail.empty());
    CHECK(q.fwd_positions.empty());

    // Exactly k bases: not too short, no skip.
    auto q2 = preprocess_query<uint16_t>(std::string(9, 'A'), 9, nullptr, config);
    CHECK_EQ(q2.skip_reason, (uint8_t)kSkipNone);
}

// Phase 1: accept_qdegen=0 with degenerate bases triggers kSkipDegenRejected.
static void test_degen_rejected() {
    std::fprintf(stderr, "-- test_degen_rejected\n");
    SearchConfig config;
    config.stage1.min_stage1_score = 1;
    config.accept_qdegen = 0;

    // Contains 'N' (degenerate)
    auto q = preprocess_query<uint16_t>("ACGTNACGTACGT", 9, nullptr, config);
    CHECK_EQ(q.skip_reason, (uint8_t)kSkipDegenRejected);

    // Pure ACGT — should not be rejected
    auto q2 = preprocess_query<uint16_t>("ACGTACGTACGTA", 9, nullptr, config);
    CHECK_EQ(q2.skip_reason, (uint8_t)kSkipNone);
}

// Phase 1: truly invalid characters (non-IUPAC, non-ACGT) trigger
// kSkipInvalidChar. '*' is not part of ncbi4na nor ATGC.
static void test_invalid_char() {
    std::fprintf(stderr, "-- test_invalid_char\n");
    SearchConfig config;
    config.stage1.min_stage1_score = 1;
    config.accept_qdegen = 1;

    auto q = preprocess_query<uint16_t>("ACGTAC*GTACGT", 9, nullptr, config);
    CHECK_EQ(q.skip_reason, (uint8_t)kSkipInvalidChar);
    // skip_detail should mention the position
    CHECK(q.skip_detail.find("pos") != std::string::npos);
}

// Phase 1: Nqkmer is the *window count*, independent of degenerate-position
// emission counts. Verify by reading back resolved_threshold derived from
// Nqkmer * P. With khx=nullptr, Nhighfreq = case 2 + case 3 only (no .khx
// exclusion = case 1 = 0). Pure ATGC has Nhighfreq=0, so threshold == ceil(Nqkmer * P).
static void test_nqkmer_window_count() {
    std::fprintf(stderr, "-- test_nqkmer_window_count\n");
    SearchConfig config;
    config.stage1.min_stage1_score = 1;
    config.min_stage1_score_frac = 0.5;
    config.strand = 2;

    // 100 bases, k=9, span=9 → Nqkmer=92, threshold = ceil(92*0.5) = 46
    auto q = preprocess_query<uint16_t>(std::string(100, 'A'), 9, nullptr, config);
    CHECK_EQ(q.skip_reason, (uint8_t)kSkipNone);
    CHECK_EQ(q.resolved_threshold_fwd, 46u);
    CHECK_EQ(q.resolved_threshold_rc, 46u);
}

// Phase 1: when threshold falls to <= 0 (impossibly low signal), the query
// is skipped with kSkipThresholdUnreachable rather than silently returning
// an empty hit list.
static void test_threshold_unreachable() {
    std::fprintf(stderr, "-- test_threshold_unreachable\n");
    SearchConfig config;
    config.stage1.min_stage1_score = 1;
    // P=1.0 is invalid for fractional mode (must be < 1.0). Use 0.99.
    config.min_stage1_score_frac = 0.99;
    config.strand = 2;
    config.max_degen_expand = 4;

    // A query mostly composed of degen-exceeded windows: strings with
    // many consecutive 'N's force expansion product to exceed max_degen_expand
    // for nearly every window position, so emitted_count ~= 0 and
    // Nhighfreq ~= Nqkmer, dropping threshold to <= 0.
    std::string seq;
    seq.append(20, 'N'); // window with 20 N's; degen product = 4^20 >> 4
    auto q = preprocess_query<uint16_t>(seq, 9, nullptr, config);
    // Either skipped (if both strands fail) or threshold==0.
    if (q.skip_reason == kSkipNone) {
        CHECK(q.resolved_threshold_fwd == 0 || q.resolved_threshold_rc == 0);
    } else {
        CHECK_EQ(q.skip_reason, (uint8_t)kSkipThresholdUnreachable);
        CHECK(q.skip_detail.find("threshold") != std::string::npos);
    }
}

// Phase 1: qlen is set on every QueryKmerData, including skip cases, so
// the OutputHit skip-marker row can carry the original query length.
static void test_qlen_populated() {
    std::fprintf(stderr, "-- test_qlen_populated\n");
    SearchConfig config;
    config.stage1.min_stage1_score = 1;
    config.accept_qdegen = 0;

    auto q1 = preprocess_query<uint16_t>("ACGTNACGTACGT", 9, nullptr, config);
    CHECK_EQ(q1.qlen, 13u);
    CHECK_EQ(q1.skip_reason, (uint8_t)kSkipDegenRejected);

    auto q2 = preprocess_query<uint16_t>(std::string(50, 'A'), 9, nullptr, config);
    CHECK_EQ(q2.qlen, 50u);
    CHECK_EQ(q2.skip_reason, (uint8_t)kSkipNone);
}

// Phase 1: skip_reason_str returns expected strings.
static void test_skip_reason_str() {
    std::fprintf(stderr, "-- test_skip_reason_str\n");
    CHECK(std::string(skip_reason_str(kSkipNone)) == "ok");
    CHECK(std::string(skip_reason_str(kSkipQueryTooShort)) == "query_too_short");
    CHECK(std::string(skip_reason_str(kSkipDegenRejected)) == "degen_rejected");
    CHECK(std::string(skip_reason_str(kSkipThresholdUnreachable)) == "threshold_unreachable");
    CHECK(std::string(skip_reason_str(kSkipInvalidChar)) == "invalid_char");
}

int main() {
    test_too_short();
    test_degen_rejected();
    test_invalid_char();
    test_nqkmer_window_count();
    test_threshold_unreachable();
    test_qlen_populated();
    test_skip_reason_str();
    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
