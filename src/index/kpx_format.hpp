#pragma once

#include <cstdint>

namespace ikafssn {

inline constexpr char KPX_MAGIC[4] = {'K', 'M', 'P', 'X'};

#pragma pack(push, 1)
struct KpxHeader {
    char     magic[4];        // 0x00: "KMPX"
    uint16_t format_version;  // 0x04
    uint8_t  k;               // 0x06
    uint8_t  reserved1;       // 0x07
    uint64_t total_postings;  // 0x08
    uint8_t  reserved2[16];   // 0x10
};
#pragma pack(pop)

static_assert(sizeof(KpxHeader) == 32, "KpxHeader must be 32 bytes");

} // namespace ikafssn
