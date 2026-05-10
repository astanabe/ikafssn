#pragma once

#include <cstdint>
#include <cstring>

namespace ikafssn {

// Elias-Fano dictionary + posting file body.  Codec selection follows
// format_version; codec_id / codec_version / block_size / tail_codec
// are reserved fields kept for header-byte stability.
inline constexpr char KIX_MAGIC[5] = {'K', 'I', 'X', '1', '1'};

// Flag bits
inline constexpr uint32_t KIX_FLAG_HAS_KSX  = 0x02; // 0=no .ksx, 1=has .ksx
inline constexpr uint32_t KIX_FLAG_OFFSET32 = 0x04; // reserved (always 0)

#pragma pack(push, 1)
struct KixHeader {
    char     magic[5];                   // 0x00: "KIX11"
    uint16_t format_version;             // 0x05
    uint8_t  k;                          // 0x07
    uint8_t  kmer_type;                  // 0x08: 0=uint16, 1=uint32
    uint32_t num_sequences;              // 0x09: number of internal sequence ids (= fragments)
    uint64_t total_distinct_postings;    // 0x0D: sum of distinct seq_id counts across all k-mers
    uint32_t flags;                      // 0x15
    uint16_t volume_index;                // 0x19
    uint16_t total_volumes;               // 0x1B
    uint16_t db_len;                      // 0x1D
    uint8_t  t;                           // 0x1F: template length (0=contiguous)
    uint8_t  template_type;               // 0x20: TemplateType enum (0=contiguous)
    char     db[32];                      // 0x21
    // ---- 32 B codec-extension area (total header size = 96 B) ----
    uint8_t  codec_id;                    // 0x41: reserved
    uint8_t  codec_version;               // 0x42: reserved
    uint16_t block_size;                  // 0x43: reserved
    uint8_t  tail_codec;                  // 0x45: reserved
    uint8_t  reserved0[3];                // 0x46
    uint32_t exception_codec_flags;       // 0x49: future use
    uint8_t  reserved1[19];               // 0x4D: pads header to a 96 B total
};
#pragma pack(pop)

static_assert(sizeof(KixHeader) == 96, "KixHeader must be 96 bytes");

} // namespace ikafssn
