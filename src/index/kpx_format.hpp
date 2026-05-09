#pragma once

#include <cstdint>

namespace ikafssn {

// .kpx v10 (Phase 1 of the fragment-indexing initiative): per-(kmer,
// seq_id) partitioned position posting file.  The dictionary and per-
// posting-list body shape are unchanged from v9 (Elias-Fano pos_offsets +
// 2-bit kind map + partition / short_occ1 / short_occ2 sub-buckets — see
// kpx_format v9 notes preserved in git history for full details).  v10
// adds the same fragment-indexing triplet to the header:
//
//   uint32_t min_seq_length    // -min_seq_length passed to ikafssnindex
//   uint32_t min_length_split  // 0 in Phase 1; non-zero from Phase 2
//   uint32_t overlap_length    // 0 in Phase 1; non-zero from Phase 2
//
// The magic widens from 4 bytes "KPX9" to 5 bytes "KPX10"; the rest of
// the header shifts +1 and one byte is taken out of the reserved tail so
// the total header size stays at 64 bytes.
inline constexpr char KPX_MAGIC[5] = {'K', 'P', 'X', '1', '0'};

// Codec selection follows format_version (since Phase 5g-1). The header
// still carries an 8-bit codec_id field for layout stability but readers
// no longer dispatch on it.

#pragma pack(push, 1)
struct KpxHeader {
    char     magic[5];                   // 0x00: "KPX10"
    uint16_t format_version;             // 0x05: 10
    uint8_t  k;                          // 0x07
    uint8_t  t;                          // 0x08: template length (0=contiguous)
    uint64_t total_position_count;       // 0x09: sum of position counts across all k-mers
    uint8_t  template_type;              // 0x11: TemplateType enum
    uint8_t  offset_type;                // 0x12: 0=uint32 offsets, 1=uint64 offsets, 0xFF=EF
    uint8_t  reserved2[14];              // 0x13
    // ---- 32 B codec-extension area (total header size = 64 B) ----
    // codec_id / codec_version / block_size / tail_codec are reserved
    // since Phase 5g-1: codec selection now follows format_version. Kept
    // as fields for header-byte stability.
    uint8_t  codec_id;                   // 0x21: 0 (reserved)
    uint8_t  codec_version;              // 0x22: 0 (reserved)
    uint16_t block_size;                 // 0x23: 0 (reserved)
    uint8_t  tail_codec;                 // 0x25: 0 (reserved)
    // v10 additions (offset 0x26, 12 bytes total)
    uint32_t min_seq_length;             // 0x26: short-sequence filter threshold
    uint32_t min_length_split;           // 0x2A: 0 in Phase 1
    uint32_t overlap_length;             // 0x2E: 0 in Phase 1
    uint8_t  reserved0[14];              // 0x32: pads header to a 64 B total
};
#pragma pack(pop)

static_assert(sizeof(KpxHeader) == 64, "KpxHeader v10 must be 64 bytes");

} // namespace ikafssn
