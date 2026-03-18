#pragma once

#include <cstdint>

namespace ikafssn {

inline constexpr char KPX_MAGIC[4] = {'K', 'M', 'P', 'X'};

#pragma pack(push, 1)
struct KpxHeader {
    char     magic[4];        // 0x00: "KMPX"
    uint16_t format_version;  // 0x04
    uint8_t  k;               // 0x06
    uint8_t  t;               // 0x07: template length (0=contiguous)
    uint64_t total_postings;  // 0x08
    uint8_t  template_type;   // 0x10: TemplateType enum value (0=contiguous)
    uint8_t  offset_type;     // 0x11: 0=uint32 offsets, 1=uint64 offsets
    uint8_t  reserved2[14];   // 0x12
};
#pragma pack(pop)

static_assert(sizeof(KpxHeader) == 32, "KpxHeader must be 32 bytes");

} // namespace ikafssn
