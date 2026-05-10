#pragma once

#include <cstdint>

namespace ikafssn {

// Bitset of k-mers excluded from the index because their distinct
// seq_id count exceeded the build-time threshold.  Magic "KMHX"
// carries no version digit; format_version is checked at open time.
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
