#include "index/fragment_splitter.hpp"

namespace ikafssn::fragment_splitter {

namespace {

// Append fragments derived from a single valid segment to `out`,
// matching the kafssstore split_long_sequence_positions() formula.
void split_segment(uint32_t seg_start_1b,
                   uint32_t seg_end_1b,
                   uint32_t min_length_split,
                   uint32_t overlap_length,
                   std::vector<Fragment>& out) {
    const uint32_t seg_len = seg_end_1b - seg_start_1b + 1;

    // nsplit = floor((seg_len - ovllen) / (minsplitlen - ovllen)).
    // The kafsss code uses Perl's `int()`, which truncates toward 0,
    // matching unsigned integer division here.  When seg_len < minsplitlen
    // (or even seg_len <= overlap_length), this floors to 0 and we emit
    // the segment as a single fragment.
    uint32_t nsplit = 0;
    if (seg_len > overlap_length && min_length_split > overlap_length) {
        nsplit = (seg_len - overlap_length) / (min_length_split - overlap_length);
    }

    if (nsplit < 1) {
        out.push_back({seg_start_1b, seg_end_1b});
        return;
    }

    // L = seg_len - (nsplit - 1) * overlap_length is the total *unique*
    // base count carved into nsplit pieces.  q / r distribute it: the
    // first r pieces get one extra base.  Each piece (except the last)
    // then gains overlap_length bases on its tail to overlap with the
    // next piece.  Cumulative offset advances by len - overlap_length
    // after each piece, leaving the trailing overlap region in place.
    const uint32_t L = seg_len - (nsplit - 1) * overlap_length;
    const uint32_t q = L / nsplit;
    const uint32_t r = L % nsplit;

    uint32_t pos = 0;  // 0-based offset within the segment
    for (uint32_t i = 0; i < nsplit; ++i) {
        uint32_t len = q + (i < r ? 1u : 0u);
        if (i + 1 != nsplit) len += overlap_length;
        const uint32_t fragment_start = seg_start_1b + pos;
        const uint32_t fragment_end   = fragment_start + len - 1;
        out.push_back({fragment_start, fragment_end});
        pos += len - overlap_length;
    }
}

} // namespace

std::vector<Fragment> split(uint32_t parent_length,
                            const std::vector<AmbiguityEntry>& ambig_entries,
                            uint32_t min_length_split,
                            uint32_t overlap_length) {
    std::vector<Fragment> fragments;

    if (parent_length == 0) {
        return fragments;
    }

    // Degenerate / splitting-disabled mode: return one fragment covering
    // the whole parent.  This intentionally ignores N runs because the
    // existing builder treats the entire parent as a single sequence.
    if (min_length_split == 0) {
        fragments.push_back({1u, parent_length});
        return fragments;
    }

    // Build the list of valid (non-N) segments by walking the ambiguity
    // table in order.  Only ncbi4na == 0xF ("any base") entries count as
    // breaks; all other codes (R, Y, W, K, M, S, B, D, H, V) are treated
    // as ambiguous-but-still-part-of-the-segment, matching kafsss DNA2
    // mode where only N is invalid.  Ambiguity entries are guaranteed
    // sorted by position by AmbiguityParser.
    uint32_t cursor_1b = 1;  // next position the segment can start from (1-based)
    for (const auto& e : ambig_entries) {
        if (e.ncbi4na != 0xF) continue;
        const uint32_t n_start_1b = e.position + 1;
        const uint32_t n_end_1b   = e.position + e.run_length;
        if (n_start_1b > cursor_1b) {
            split_segment(cursor_1b, n_start_1b - 1,
                          min_length_split, overlap_length, fragments);
        }
        if (n_end_1b + 1 > cursor_1b) {
            cursor_1b = n_end_1b + 1;
        }
    }
    if (cursor_1b <= parent_length) {
        split_segment(cursor_1b, parent_length,
                      min_length_split, overlap_length, fragments);
    }

    // If the parent is entirely N (or empty after filtering), fall back
    // to a single fragment covering the parent.  This keeps the SeqId
    // accounting trivial (every parent contributes >= 1 fragment) and
    // matches the degenerate / "no valid bases" fallback in kafsss
    // (which still records the OID).
    if (fragments.empty()) {
        fragments.push_back({1u, parent_length});
    }

    return fragments;
}

} // namespace ikafssn::fragment_splitter
