#pragma once

#include <cstdint>

namespace ikafssn {

// .kpx v5 (Phase 5g-1): position posting format with custom FOR-within-block.
// Stores absolute positions; each 128-element block subtracts its min before
// bitpacking. Data layout is unchanged from v4; magic and format_version
// bump with the rest of the index family for alignment (Phase 5c policy).
inline constexpr char KPX_MAGIC[4] = {'K', 'P', 'X', '5'};

// Codec selection moved to format_version. The header still carries an
// 8-bit codec_id field for layout stability but readers no longer dispatch
// on it.

#pragma pack(push, 1)
struct KpxHeader {
    char     magic[4];                   // 0x00: "KPX5"
    uint16_t format_version;             // 0x04: 5
    uint8_t  k;                          // 0x06
    uint8_t  t;                          // 0x07: template length (0=contiguous)
    uint64_t total_postings;             // 0x08
    uint8_t  template_type;              // 0x10: TemplateType enum
    uint8_t  offset_type;                // 0x11: 0=uint32 offsets, 1=uint64 offsets
    uint8_t  reserved2[14];              // 0x12
    // ---- 32 B codec-extension area (total header size = 64 B) ----
    // codec_id / codec_version / block_size / tail_codec are reserved
    // since Phase 5g-1: codec selection now follows format_version. Kept
    // as fields for header-byte stability.
    uint8_t  codec_id;                   // 0x20: 0 (reserved)
    uint8_t  codec_version;              // 0x21: 0 (reserved)
    uint16_t block_size;                 // 0x22: 0 (reserved)
    uint8_t  tail_codec;                 // 0x24: 0 (reserved)
    uint8_t  reserved0[27];              // 0x25
};
#pragma pack(pop)

static_assert(sizeof(KpxHeader) == 64, "KpxHeader v5 must be 64 bytes");

} // namespace ikafssn
