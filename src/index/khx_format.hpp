#pragma once

#include <cstdint>

namespace ikafssn {

// .khx v9 (Phase 7): data layout is unchanged from v3; the
// format_version field bumps along with the rest of the index family
// for alignment (Phase 5c policy).  Magic stays "KMHX" — see the
// comment in src/index/ksx_format.hpp for why .ksx / .khx do not
// carry a version digit in the magic string.  The semantic of the
// bitset is "k-mer i was excluded because its **distinct seq_id
// count** exceeded the build-time threshold" (matching the .kix
// distinct-seq_id semantic established in v7).
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
