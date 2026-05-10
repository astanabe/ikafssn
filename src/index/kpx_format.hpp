#pragma once

#include <cstdint>

namespace ikafssn {

// .kpx v10: per-(kmer, seq_id) partitioned position posting file. The
// dictionary uses Elias-Fano pos_offsets, and each posting list body is
// laid out as a 2-bit kind map + partition / short_occ1 / short_occ2
// sub-buckets. Header carries the fragment-indexing triplet (min_seq_length,
// min_length_split, overlap_length) plus the standard codec-extension area.
// Codec selection follows format_version; codec_id / codec_version /
// block_size / tail_codec are reserved fields kept for header-byte stability.
inline constexpr char KPX_MAGIC[5] = {'K', 'P', 'X', '1', '0'};

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
    uint8_t  codec_id;                   // 0x21: reserved
    uint8_t  codec_version;              // 0x22: reserved
    uint16_t block_size;                 // 0x23: reserved
    uint8_t  tail_codec;                 // 0x25: reserved
    uint32_t min_seq_length;             // 0x26: short-sequence filter threshold
    uint32_t min_length_split;           // 0x2A: fragment splitter threshold (0 = no split)
    uint32_t overlap_length;             // 0x2E: overlap between adjacent fragments
    uint8_t  reserved0[14];              // 0x32: pads header to a 64 B total
};
#pragma pack(pop)

static_assert(sizeof(KpxHeader) == 64, "KpxHeader v10 must be 64 bytes");

} // namespace ikafssn
