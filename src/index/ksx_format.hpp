#pragma once

#include <cstdint>

namespace ikafssn {

// .ksx v10 (Phase 1 of the fragment-indexing initiative): the .ksx is now
// a two-stage table.  The first stage records each parent OID's full-
// length metadata (parent length, BLAST DB volume-local OID, accession);
// the second stage records the fragments derived from each parent.  In
// Phase 1 each parent has exactly one fragment spanning the whole parent
// (`fragment_start = 1`, `fragment_end = parent_length`), so num_parents
// equals num_sequences.  Phase 2 introduces real splitting.
//
// On-disk layout (after this 32 B header):
//   parent_lengths[num_parents]            : u32
//   parent_blast_oids[num_parents]         : u32
//   parent_acc_offsets[num_parents + 1]    : u32
//   parent_acc_strings                     : concatenated bytes (no NUL)
//   fragment_parent_idx[num_sequences]     : u32
//   fragment_start[num_sequences]          : u32 (1-based, parent-relative)
//   fragment_end[num_sequences]            : u32 (1-based, parent-relative)
//
// Magic stays "KMSX" — the version-digit-suffix convention is reserved
// for the codec-bearing .kix and .kpx files; .ksx / .khx use only
// format_version for version detection.
inline constexpr char KSX_MAGIC[4] = {'K', 'M', 'S', 'X'};

#pragma pack(push, 1)
struct KsxHeader {
    char     magic[4];          // 0x00: "KMSX"
    uint16_t format_version;    // 0x04: 10
    uint16_t reserved1;         // 0x06
    uint32_t num_sequences;     // 0x08: number of internal sequence ids (= fragments)
    // v10 additions
    uint32_t min_seq_length;    // 0x0C: short-sequence filter threshold
    uint32_t min_length_split;  // 0x10: 0 in Phase 1
    uint32_t overlap_length;    // 0x14: 0 in Phase 1
    uint32_t num_parents;       // 0x18: number of parent OIDs
    uint8_t  reserved2[4];      // 0x1C: pads header to 32 B total
};
#pragma pack(pop)

static_assert(sizeof(KsxHeader) == 32, "KsxHeader v10 must be 32 bytes");

} // namespace ikafssn
