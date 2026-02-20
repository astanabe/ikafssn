#pragma once

#include <cstdint>

namespace ikafssn {

inline constexpr char KSX_MAGIC[4] = {'K', 'M', 'S', 'X'};

#pragma pack(push, 1)
struct KsxHeader {
    char     magic[4];        // 0x00: "KMSX"
    uint16_t format_version;  // 0x04
    uint16_t reserved1;       // 0x06
    uint32_t num_sequences;   // 0x08
    uint8_t  reserved2[20];   // 0x0C
};
#pragma pack(pop)

static_assert(sizeof(KsxHeader) == 32, "KsxHeader must be 32 bytes");

} // namespace ikafssn
