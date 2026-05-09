// Phase 5 unit test: confirm src/index/fragment_splitter::split() reproduces
// kafssstore.pl's split_long_sequence_positions (DNA2 mode, ncbi4na==0xF
// cuts, calcsegment2 formula).  We anchor the C++ port both to fixed
// expected outputs (so a regression in either the segment scan or the
// per-segment split would fail loudly) and to a stand-alone Perl-equivalent
// reimplementation in C++ (so larger random parents are also covered).
//
// The reference implementation lives at
//   ~/claudecode/kafsss/kafssstore.pl:1078-1180
//
// Coverage notes on the kafsss formula:
//   * nsplit = floor((seg_len - ovllen) / (minsplitlen - ovllen)).
//   * If nsplit < 1 the segment becomes one fragment unchanged.
//   * Else L = seg_len - (nsplit - 1) * ovllen, q = L / nsplit,
//     r = L % nsplit; the leftmost r fragments get one extra base.
//   * Each fragment except the last gains ovllen trailing bases that
//     overlap the next fragment.  As a deliberate consequence of this
//     formula, the LAST fragment ends at L (not at seg_len) when
//     ovllen > 0 and nsplit > 1, leaving (nsplit-1)*ovllen bases at
//     the parent tail uncovered.  This is the kafsss spec; the C++
//     port preserves it byte-for-byte.

#include "test_util.hpp"
#include "core/ambiguity_parser.hpp"
#include "index/fragment_splitter.hpp"

#include <cstdint>
#include <random>
#include <string>
#include <vector>

using namespace ikafssn;
using namespace ikafssn::fragment_splitter;

namespace {

// Stand-alone re-implementation of split_long_sequence_positions for
// random / property-based comparison.  Walks `valid_chars`-style text
// to identify segments, then carves each segment with the calcsegment2
// formula used in the Perl reference.  We model "N runs" by passing
// segment boundaries directly so we can drive both the C++ port and the
// reference from the same parent string.
struct PerlFragment {
    uint32_t start;
    uint32_t end;
};

std::vector<PerlFragment> reference_split_segments(
    const std::vector<std::pair<uint32_t, uint32_t>>& segments,
    uint32_t min_length_split, uint32_t overlap_length) {
    std::vector<PerlFragment> out;
    for (auto [seg_start, seg_end] : segments) {
        const uint32_t seg_len = seg_end - seg_start + 1;
        uint32_t nsplit = 0;
        if (seg_len > overlap_length && min_length_split > overlap_length) {
            nsplit = (seg_len - overlap_length) / (min_length_split - overlap_length);
        }
        if (nsplit < 1) {
            out.push_back({seg_start, seg_end});
            continue;
        }
        const uint32_t L = seg_len - (nsplit - 1) * overlap_length;
        const uint32_t q = L / nsplit;
        const uint32_t r = L % nsplit;
        uint32_t pos = 0;
        for (uint32_t i = 0; i < nsplit; ++i) {
            uint32_t len = q + (i < r ? 1u : 0u);
            if (i + 1 != nsplit) len += overlap_length;
            const uint32_t fragment_start = seg_start + pos;
            const uint32_t fragment_end   = fragment_start + len - 1;
            out.push_back({fragment_start, fragment_end});
            pos += len - overlap_length;
        }
    }
    return out;
}

// Build an AmbiguityEntry list that cuts `parent_length` into the given
// valid `segments` (1-based, inclusive).  Each gap between consecutive
// segments becomes one ncbi4na==0xF entry.  Trailing N runs after the
// last segment are also emitted so the splitter sees the full ambiguity
// table.
std::vector<AmbiguityEntry> ambig_for_segments(
    uint32_t parent_length,
    const std::vector<std::pair<uint32_t, uint32_t>>& segments) {
    std::vector<AmbiguityEntry> out;
    uint32_t cursor = 1;
    for (auto [s, e] : segments) {
        if (s > cursor) {
            AmbiguityEntry a;
            a.position   = cursor - 1;        // 0-based start of N run
            a.run_length = s - cursor;        // number of N bases
            a.ncbi4na    = 0xF;
            out.push_back(a);
        }
        cursor = e + 1;
    }
    if (cursor <= parent_length) {
        AmbiguityEntry a;
        a.position   = cursor - 1;
        a.run_length = parent_length - cursor + 1;
        a.ncbi4na    = 0xF;
        out.push_back(a);
    }
    return out;
}

void check_fragments_eq(const std::vector<Fragment>& got,
                        const std::vector<PerlFragment>& expected,
                        const char* label) {
    CHECK_EQ(got.size(), expected.size());
    if (got.size() != expected.size()) {
        std::fprintf(stderr, "  [%s] size mismatch\n", label);
        return;
    }
    for (size_t i = 0; i < got.size(); ++i) {
        if (got[i].start != expected[i].start ||
            got[i].end   != expected[i].end) {
            std::fprintf(stderr,
                "  [%s] fragment %zu: got [%u..%u], expected [%u..%u]\n",
                label, i, got[i].start, got[i].end,
                expected[i].start, expected[i].end);
            g_fail_count++;
        } else {
            g_pass_count++;
        }
    }
}

}  // namespace

// 1. Degenerate / splitting-disabled: min_length_split == 0 always emits one
// fragment spanning the whole parent, regardless of ambiguity entries.
static void test_degenerate_disables_splitting() {
    std::fprintf(stderr, "-- test_degenerate_disables_splitting\n");

    // A parent with mid-length N runs that would otherwise cut into pieces.
    auto ambig = ambig_for_segments(1000, {{1, 400}, {501, 1000}});
    auto frags = split(1000, ambig, /*min_length_split*/ 0, /*overlap*/ 0);
    CHECK_EQ(frags.size(), 1u);
    CHECK_EQ(frags[0].start, 1u);
    CHECK_EQ(frags[0].end,   1000u);

    // Even with overlap_length set, min_length_split==0 still means
    // "splitting disabled" by the public contract.
    frags = split(1000, ambig, 0, 200);
    CHECK_EQ(frags.size(), 1u);
}

// 2. Short parent (no segments after N filtering): fall back to a single
// fragment covering the entire parent.  Mirrors the kafsss "no valid
// bases" fallback (the OID is still recorded).
static void test_all_n_parent() {
    std::fprintf(stderr, "-- test_all_n_parent\n");
    AmbiguityEntry a{0u, 500u, 0xF};
    auto frags = split(500, {a}, /*min_length_split*/ 1500, /*overlap*/ 200);
    CHECK_EQ(frags.size(), 1u);
    CHECK_EQ(frags[0].start, 1u);
    CHECK_EQ(frags[0].end,   500u);
}

// 3. Non-N ambiguity codes (R/Y/W/...) MUST stay inside the segment.
// kafsss DNA2 mode only treats ncbi4na==0xF (N) as a segment break.
static void test_non_N_codes_do_not_break_segment() {
    std::fprintf(stderr, "-- test_non_N_codes_do_not_break_segment\n");
    // Parent of 5000 bp with several non-N degenerate codes scattered.
    std::vector<AmbiguityEntry> ambig = {
        {  100u, 1u, 0x5 /* R */},
        {  500u, 1u, 0xA /* Y */},
        { 2000u, 4u, 0x9 /* W run */},
        { 4500u, 1u, 0x6 /* M */},
    };
    auto frags = split(5000, ambig, 1500, 200);

    // Reference: one segment [1, 5000], then split with the calcsegment2
    // formula.  nsplit = (5000 - 200) / (1500 - 200) = 4800 / 1300 = 3.
    auto expected = reference_split_segments({{1, 5000}}, 1500, 200);
    CHECK_EQ(expected.size(), 3u);
    check_fragments_eq(frags, expected, "non-N codes");
}

// 4. Single N run cuts the parent into two valid segments; each is then
// independently split.  This is the canonical kafsss DNA2 case.
static void test_one_n_run_cuts_into_two_segments() {
    std::fprintf(stderr, "-- test_one_n_run_cuts_into_two_segments\n");
    // Segments: [1, 2000] and [2010, 6000]; the N run is positions 2001-2009.
    AmbiguityEntry a{2000u, 9u, 0xF};
    auto frags = split(6000, {a}, 1500, 200);

    auto expected = reference_split_segments(
        {{1, 2000}, {2010, 6000}}, 1500, 200);
    check_fragments_eq(frags, expected, "one N run");
}

// 5. nsplit < 1 (segment shorter than min_length_split): emit segment
// untouched even if it falls under min_length_split.
static void test_segment_shorter_than_min_split() {
    std::fprintf(stderr, "-- test_segment_shorter_than_min_split\n");
    // Parent of 1000 bp, no N run.  min_length_split = 1500 -> nsplit = 0.
    auto frags = split(1000, {}, 1500, 200);
    CHECK_EQ(frags.size(), 1u);
    CHECK_EQ(frags[0].start, 1u);
    CHECK_EQ(frags[0].end,   1000u);
}

// 6. Exact LONGCHR_1 / LONGCHR_2 expected output for the SSU_longdb
// fixture parameters used by the integration tests below.  Hard-coded
// numbers come straight from running the Perl reference by hand.
static void test_longchr_expected_split() {
    std::fprintf(stderr, "-- test_longchr_expected_split\n");

    // LONGCHR_1: 5171 bp, no N runs (concatenated SSU sequences).
    // nsplit = (5171 - 200) / (1500 - 200) = 4971 / 1300 = 3.
    // L = 5171 - 2*200 = 4771; q = 1590; r = 1.
    // i=0: len = 1591 + 200 = 1791  -> [1, 1791]
    // i=1: len = 1590 + 200 = 1790  -> [1592, 3381]
    // i=2: len = 1590                -> [3182, 4771]
    auto frags = split(5171u, {}, 1500u, 200u);
    CHECK_EQ(frags.size(), 3u);
    if (frags.size() == 3) {
        CHECK_EQ(frags[0].start, 1u);     CHECK_EQ(frags[0].end, 1791u);
        CHECK_EQ(frags[1].start, 1592u);  CHECK_EQ(frags[1].end, 3381u);
        CHECK_EQ(frags[2].start, 3182u);  CHECK_EQ(frags[2].end, 4771u);
    }

    // LONGCHR_2: 5195 bp.  nsplit = 4995 / 1300 = 3.
    // L = 5195 - 400 = 4795; q = 1598; r = 1.
    // i=0: len = 1599 + 200 = 1799  -> [1, 1799]
    // i=1: len = 1598 + 200 = 1798  -> [1600, 3397]
    // i=2: len = 1598                -> [3198, 4795]
    frags = split(5195u, {}, 1500u, 200u);
    CHECK_EQ(frags.size(), 3u);
    if (frags.size() == 3) {
        CHECK_EQ(frags[0].start, 1u);     CHECK_EQ(frags[0].end, 1799u);
        CHECK_EQ(frags[1].start, 1600u);  CHECK_EQ(frags[1].end, 3397u);
        CHECK_EQ(frags[2].start, 3198u);  CHECK_EQ(frags[2].end, 4795u);
    }
}

// 7. Property test: random parent lengths and N runs, compared against the
// stand-alone reference re-implementation.  Catches divergence between
// the C++ port and the Perl formula in corner regimes that the fixed
// cases don't exercise (e.g., exact-multiple seg_len, segment of length
// overlap+1, multiple back-to-back N runs).
static void test_property_random_against_reference() {
    std::fprintf(stderr, "-- test_property_random_against_reference\n");
    std::mt19937 rng(20260510u);
    std::uniform_int_distribution<uint32_t> parent_dist(100u, 12000u);
    std::uniform_int_distribution<uint32_t> minsplit_dist(500u, 2500u);
    std::uniform_int_distribution<uint32_t> overlap_frac(2u, 8u);   // 1/8..1/2
    std::uniform_int_distribution<uint32_t> n_run_count(0u, 4u);
    std::uniform_int_distribution<uint32_t> n_run_len(1u, 200u);

    int trials = 100;
    for (int t = 0; t < trials; ++t) {
        const uint32_t parent_length = parent_dist(rng);
        const uint32_t minsplit = minsplit_dist(rng);
        // overlap ranges [1, minsplit/2 - 1] so the splitter precondition holds.
        const uint32_t denom = overlap_frac(rng);
        uint32_t overlap = minsplit / denom;
        if (overlap == 0) overlap = 1;
        if (overlap >= minsplit) overlap = minsplit - 1;
        if (overlap > minsplit / 2) overlap = minsplit / 2;
        if (overlap == 0) overlap = 1;

        // Sprinkle a few N runs at random non-overlapping positions.
        const uint32_t nruns = n_run_count(rng);
        std::vector<AmbiguityEntry> ambig;
        std::vector<std::pair<uint32_t, uint32_t>> n_intervals;
        for (uint32_t i = 0; i < nruns; ++i) {
            uint32_t pos = std::uniform_int_distribution<uint32_t>(
                0, parent_length > 1 ? parent_length - 1 : 0)(rng);
            uint32_t len = n_run_len(rng);
            if (pos + len > parent_length) len = parent_length - pos;
            if (len == 0) continue;
            n_intervals.push_back({pos, pos + len - 1});  // 0-based inclusive
        }
        // Sort and merge overlapping N runs so the segment derivation below
        // mirrors what the splitter sees.  AmbiguityParser already sorts by
        // position; we additionally collapse overlaps for correctness.
        std::sort(n_intervals.begin(), n_intervals.end());
        std::vector<std::pair<uint32_t, uint32_t>> merged;
        for (auto& iv : n_intervals) {
            if (!merged.empty() && iv.first <= merged.back().second + 1) {
                merged.back().second = std::max(merged.back().second, iv.second);
            } else {
                merged.push_back(iv);
            }
        }
        for (auto& iv : merged) {
            AmbiguityEntry a;
            a.position   = iv.first;
            a.run_length = iv.second - iv.first + 1;
            a.ncbi4na    = 0xF;
            ambig.push_back(a);
        }

        // Derive valid (non-N) segments in 1-based, inclusive coords.
        std::vector<std::pair<uint32_t, uint32_t>> segments;
        uint32_t cursor = 1;
        for (auto& iv : merged) {
            const uint32_t n_start_1b = iv.first + 1;
            const uint32_t n_end_1b   = iv.second + 1;
            if (n_start_1b > cursor) segments.push_back({cursor, n_start_1b - 1});
            cursor = n_end_1b + 1;
        }
        if (cursor <= parent_length) segments.push_back({cursor, parent_length});

        auto got = split(parent_length, ambig, minsplit, overlap);
        std::vector<PerlFragment> expected;
        if (segments.empty()) {
            expected.push_back({1u, parent_length});
        } else {
            expected = reference_split_segments(segments, minsplit, overlap);
        }
        char label[64];
        std::snprintf(label, sizeof(label),
                      "trial %d (plen=%u min=%u ovl=%u nrun=%zu)",
                      t, parent_length, minsplit, overlap, merged.size());
        check_fragments_eq(got, expected, label);
    }
}

int main() {
    test_degenerate_disables_splitting();
    test_all_n_parent();
    test_non_N_codes_do_not_break_segment();
    test_one_n_run_cuts_into_two_segments();
    test_segment_shorter_than_min_split();
    test_longchr_expected_split();
    test_property_random_against_reference();
    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
