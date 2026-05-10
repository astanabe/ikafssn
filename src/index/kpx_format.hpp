#pragma once

#include <cstdint>

namespace ikafssn {

// Per-(kmer, seq_id) partitioned position posting file.  The
// dictionary uses Elias-Fano pos_offsets; each posting list body is
// laid out as a 2-bit kind map plus partition / short_occ1 /
// short_occ2 sub-buckets.  Codec selection follows format_version;
// codec_id / codec_version / block_size / tail_codec are reserved
// fields kept for header-byte stability.
inline constexpr char KPX_MAGIC[5] = {'K', 'P', 'X', '1', '1'};

#pragma pack(push, 1)
struct KpxHeader {
    char     magic[5];                   // 0x00: "KPX11"
    uint16_t format_version;             // 0x05
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
    uint8_t  reserved0[26];              // 0x26: pads header to a 64 B total
};
#pragma pack(pop)

static_assert(sizeof(KpxHeader) == 64, "KpxHeader must be 64 bytes");

} // namespace ikafssn
