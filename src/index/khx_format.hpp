#pragma once

#include <cstdint>

namespace ikafssn {

inline constexpr char KHX_MAGIC[4] = {'K', 'M', 'H', 'X'};

#pragma pack(push, 1)
struct KhxHeader {
    char     magic[4];        // 0x00: "KMHX"
    uint16_t format_version;  // 0x04
    uint8_t  k;               // 0x06
    uint8_t  t;               // 0x07: template length (0=contiguous)
    uint8_t  template_type;   // 0x08: TemplateType enum value (0=contiguous)
    uint8_t  reserved2[23];   // 0x09
};
#pragma pack(pop)

static_assert(sizeof(KhxHeader) == 32, "KhxHeader must be 32 bytes");

} // namespace ikafssn
