#pragma once

#include <cstdint>

namespace ikafssn {

// Two-stage parent / fragment metadata table.  The parent stage
// records each parent OID's `(parent_length, blast_oid, accession)`;
// the fragment stage records each fragment as `(parent_idx,
// fragment_start, fragment_end)`.  When fragment splitting is disabled
// every parent has exactly one fragment spanning the whole parent.
//
// On-disk layout (after the 24 B header):
//   parent_lengths[num_parents]            : u32
//   parent_blast_oids[num_parents]         : u32
//   parent_acc_offsets[num_parents + 1]    : u32
//   parent_acc_strings                     : concatenated bytes (no NUL)
//   fragment_parent_idx[num_sequences]     : u32
//   fragment_start[num_sequences]          : u32 (1-based, parent-relative)
//   fragment_end[num_sequences]            : u32 (1-based, parent-relative)
//
// Magic "KMSX" carries no version digit; .ksx / .khx rely on
// format_version alone for version detection.
inline constexpr char KSX_MAGIC[4] = {'K', 'M', 'S', 'X'};

#pragma pack(push, 1)
struct KsxHeader {
    char     magic[4];          // 0x00: "KMSX"
    uint16_t format_version;    // 0x04
    uint16_t reserved1;         // 0x06
    uint32_t num_sequences;     // 0x08: number of internal sequence ids (= fragments)
    uint32_t num_parents;       // 0x0C: number of parent OIDs
    uint8_t  reserved2[8];      // 0x10: pads header to 24 B total
};
#pragma pack(pop)

static_assert(sizeof(KsxHeader) == 24, "KsxHeader must be 24 bytes");

} // namespace ikafssn
