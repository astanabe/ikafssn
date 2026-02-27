#pragma once

#include <cstdint>
#include <cstring>

namespace ikafssn {

inline constexpr char KIX_MAGIC[4] = {'K', 'M', 'I', 'X'};

// Flag bits
inline constexpr uint32_t KIX_FLAG_SEQ_ID_WIDTH = 0x01; // 0=uint32, 1=uint64 (future)
inline constexpr uint32_t KIX_FLAG_HAS_KSX      = 0x02; // 0=no .ksx, 1=has .ksx

#pragma pack(push, 1)
struct KixHeader {
    char     magic[4];        // 0x00: "KMIX"
    uint16_t format_version;  // 0x04
    uint8_t  k;               // 0x06
    uint8_t  kmer_type;       // 0x07: 0=uint16, 1=uint32
    uint32_t num_sequences;   // 0x08
    uint64_t total_postings;  // 0x0C
    uint32_t flags;           // 0x14
    uint16_t volume_index;    // 0x18
    uint16_t total_volumes;   // 0x1A
    uint16_t db_len;          // 0x1C
    uint8_t  reserved[2];     // 0x1E
    char     db[32];          // 0x20
};
#pragma pack(pop)

static_assert(sizeof(KixHeader) == 64, "KixHeader must be 64 bytes");

} // namespace ikafssn
