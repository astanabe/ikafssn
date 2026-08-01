#pragma once

// Conversion between ikafssn's internal hit coordinates and the coordinate
// convention BLAST's `-outfmt 6` uses, which is what the TSV / JSON writers
// emit and the TSV reader consumes.
//
// Internal (OutputHit / ResponseHit / ChainResult):
//   - 0-based half-open, [start, end)
//   - query positions are query-strand-relative: on the reverse strand they
//     index the reverse complement of the query, not the query itself
//   - subject positions are parent-relative and always ascending
//
// Output (BLAST `-outfmt 6`):
//   - 1-based inclusive, [start, end]
//   - query positions are query-relative, so qstart < qend on both strands
//   - subject positions are parent-relative in forward numbering but run in
//     alignment order, so a reverse-strand hit has sstart > send
//
// SAM/BAM output does not go through here: htslib's POS is 0-based and takes
// the internal sstart directly.

#include <cstdint>

namespace ikafssn {

struct OutputCoords {
    uint32_t qstart;
    uint32_t qend;
    uint32_t sstart;
    uint32_t send;
};

inline OutputCoords to_output_coords(uint32_t q_start, uint32_t q_end,
                                     uint32_t s_start, uint32_t s_end,
                                     uint32_t qlen, bool is_reverse) {
    if (!is_reverse) {
        return OutputCoords{q_start + 1, q_end, s_start + 1, s_end};
    }
    // Fold the query back into the forward frame, which reverses the interval,
    // and report the subject in alignment order (descending).
    return OutputCoords{qlen - q_end + 1, qlen - q_start, s_end, s_start + 1};
}

inline void to_internal_coords(uint32_t qstart, uint32_t qend,
                               uint32_t sstart, uint32_t send,
                               uint32_t qlen, bool is_reverse,
                               uint32_t& q_start, uint32_t& q_end,
                               uint32_t& s_start, uint32_t& s_end) {
    if (!is_reverse) {
        q_start = qstart - 1;
        q_end   = qend;
        s_start = sstart - 1;
        s_end   = send;
    } else {
        q_start = qlen - qend;
        q_end   = qlen - qstart + 1;
        s_start = send - 1;
        s_end   = sstart;
    }
}

} // namespace ikafssn
