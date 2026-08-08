#pragma once

#include <cstdint>
#include <vector>

#include "core/ambiguity_parser.hpp"

namespace ikafssn::fragment_splitter {

// One overlapping fragment derived from a parent OID.  Coordinates are
// 1-based, inclusive, and parent-relative.  Adjacent fragments belonging
// to the same valid (non-N) segment overlap by `overlap_length` bases.
struct Fragment {
    uint32_t start;   // 1-based, inclusive
    uint32_t end;     // 1-based, inclusive (start <= end <= parent_length)
};

// Port of kafssstore's split_long_sequence_positions() (DNA2 mode).
//
// Given a parent of length `parent_length` and its BLAST-DB ambiguity
// table, split the parent into overlapping fragments of approximately
// `min_length_split` bases (with `overlap_length` shared between
// adjacent fragments of the same valid segment).
//
// Splitting policy:
//
//   * Ambiguity entries with ncbi4na == 0xF (i.e. "N" runs) cut the
//     parent into "valid segments".  Other ambiguity codes
//     (R, Y, W, ...) do *not* break a segment — they remain inside it
//     and are handled downstream by the k-mer scanner's degenerate
//     expansion path.
//
//   * Each valid segment is split using the kafsss/calcsegment2 formula:
//
//         nsplit = (seg_len - overlap_length) / (min_length_split - overlap_length)
//
//     If `nsplit < 1`, the segment is emitted as a single fragment
//     (even if it's shorter than `min_length_split`).  Otherwise the
//     segment is divided into `nsplit` fragments whose lengths sum to
//     `seg_len + (nsplit - 1) * overlap_length`, rounded so the
//     leftmost `r` fragments are one base longer.
//
//   * `min_length_split == 0` is a degenerate / "splitting disabled"
//     value: the function returns a single fragment spanning
//     `[1, parent_length]`, regardless of N runs (callers that want
//     N-aware degenerate behaviour should special-case this themselves).
//
// `overlap_length` is required to be < min_length_split / 2 when
// splitting is enabled; the caller is responsible for that
// precondition (the CLI parser in ikafssnindex enforces it).
//
// Returns at least one fragment (the degenerate single-fragment case
// applies for empty-or-all-N parents too: the caller should ensure
// `parent_length >= min_seq_length` before reaching this code).
std::vector<Fragment> split(
    uint32_t parent_length,
    const std::vector<AmbiguityEntry>& ambig_entries,
    uint32_t min_length_split,
    uint32_t overlap_length);

} // namespace ikafssn::fragment_splitter
