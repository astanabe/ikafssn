#pragma once

// Conversion between the internal hit coordinates and the BLAST `-outfmt 6`
// convention, applied by the TSV / JSON writers and the TSV reader.
//
// Internal (OutputHit / ResponseHit / ChainResult) is 0-based half-open.
// Query positions are query-strand-relative: on the reverse strand they index
// the reverse complement of the query.  Subject positions are parent-relative
// and ascending.
//
// Output is 1-based inclusive.  Query positions are query-relative, so
// qstart < qend on both strands.  Subject positions keep forward numbering but
// run in alignment order, so a reverse-strand hit has sstart > send.
//
// SAM/BAM does not use this: htslib's POS is 0-based and takes sstart as is.

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
    // Fold the query into the forward frame and descend the subject.
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
