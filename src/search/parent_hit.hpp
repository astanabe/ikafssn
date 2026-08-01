#pragma once

// Parent-side view of a Stage 2 chain.
//
// Chains are produced against fragments, in fragment-relative coordinates,
// while every consumer downstream — the output writers, Stage 3, the dedup —
// reports the parent OID instead: its accession, its length, its BLAST DB OID,
// and coordinates shifted by (fragment_start - 1).

#include <cstdint>
#include <string_view>

#include "core/types.hpp"
#include "index/ksx_reader.hpp"

namespace ikafssn {

struct ParentHit {
    std::string_view sseqid;  // parent accession
    uint32_t sstart;          // parent-relative chain start
    uint32_t send;            // parent-relative chain end
    uint32_t slen;            // parent length
    uint32_t oid;             // parent BLAST DB OID
};

// Accession and length of the parent a fragment belongs to.
struct ParentInfo {
    std::string_view sseqid;
    uint32_t slen;
};

// Shift from fragment-relative to parent-relative coordinates.
inline uint32_t fragment_shift(const KsxReader& ksx, SeqId seq_id) {
    return ksx.fragment_start(seq_id) - 1u;
}

inline ParentInfo resolve_parent_info(const KsxReader& ksx, SeqId seq_id) {
    const uint32_t parent_idx = ksx.parent_index(seq_id);
    return ParentInfo{ksx.parent_accession(parent_idx),
                      ksx.parent_length(parent_idx)};
}

inline ParentHit resolve_parent_hit(const KsxReader& ksx, const ChainResult& cr) {
    const uint32_t parent_idx = ksx.parent_index(cr.seq_id);
    const uint32_t shift = fragment_shift(ksx, cr.seq_id);
    return ParentHit{
        ksx.parent_accession(parent_idx),
        cr.s_start + shift,
        cr.s_end + shift,
        ksx.parent_length(parent_idx),
        ksx.blast_oid(parent_idx),
    };
}

}  // namespace ikafssn
