#pragma once

#include <cstdint>
#include <cstring>

namespace ikafssn {

// .kix v4: SIMD-FastPFOR* + Simple-8b posting format (Phase 5b).
inline constexpr char KIX_MAGIC[4] = {'K', 'I', 'X', '4'};

// Flag bits
inline constexpr uint32_t KIX_FLAG_SEQ_ID_WIDTH = 0x01; // 0=uint32, 1=uint64 (future)
inline constexpr uint32_t KIX_FLAG_HAS_KSX      = 0x02; // 0=no .ksx, 1=has .ksx
inline constexpr uint32_t KIX_FLAG_OFFSET32     = 0x04; // 0=uint64 offsets, 1=uint32 offsets

// Codec id (codec_id field in extended v4 header).
inline constexpr uint8_t  KIX_CODEC_PFOR_S8B    = 1; // SIMD-FastPFOR* + Simple-8b/VByte tail
// Tail codec id (tail_codec field).
inline constexpr uint8_t  KIX_TAIL_VBYTE        = 1; // FastPFor CompositeCodec inner tail (VariableByte)

#pragma pack(push, 1)
struct KixHeader {
    char     magic[4];                   // 0x00: "KIX4"
    uint16_t format_version;             // 0x04: 4
    uint8_t  k;                          // 0x06
    uint8_t  kmer_type;                  // 0x07: 0=uint16, 1=uint32
    uint32_t num_sequences;              // 0x08
    uint64_t total_postings;             // 0x0C
    uint32_t flags;                      // 0x14
    uint16_t volume_index;               // 0x18
    uint16_t total_volumes;              // 0x1A
    uint16_t db_len;                     // 0x1C
    uint8_t  t;                          // 0x1E: template length (0=contiguous)
    uint8_t  template_type;              // 0x1F: TemplateType enum (0=contiguous)
    char     db[32];                     // 0x20
    // ---- v4 codec extension (32 B, total header size = 96 B) ----
    uint8_t  codec_id;                   // 0x40: KIX_CODEC_PFOR_S8B
    uint8_t  codec_version;              // 0x41: 1
    uint16_t block_size;                 // 0x42: 128
    uint8_t  tail_codec;                 // 0x44: KIX_TAIL_VBYTE
    uint8_t  reserved0[3];               // 0x45
    uint32_t exception_codec_flags;      // 0x48: future use
    uint8_t  reserved1[20];              // 0x4C: pads header to a 96 B total
};
#pragma pack(pop)

static_assert(sizeof(KixHeader) == 96, "KixHeader v4 must be 96 bytes");

} // namespace ikafssn
