#pragma once

#include <cstdint>

namespace ikafssn {

// .kpx v4: SIMD-FastPFOR* + Simple-8b position posting format (Phase 5b).
// Stores absolute positions (no within-seq delta reset); FastPFor's per-block
// bit-width adapts to the position range so dense per-sequence clusters
// compress efficiently.
inline constexpr char KPX_MAGIC[4] = {'K', 'P', 'X', '4'};

// codec_id values (codec_id field).
//   1 — Phase 5b: FastPFor CompositeCodec on absolute positions (legacy).
//   2 — Phase 5e: custom FOR-within-block — each 128-element block
//       subtracts its min before bitpacking, dramatically lowering the
//       effective bit-width on long-sequence DBs (e.g. NCBI nt-class).
//       Tail uses a single FOR base + varint stream.
inline constexpr uint8_t KPX_CODEC_PFOR_S8B = 1; // legacy
inline constexpr uint8_t KPX_CODEC_PFOR_FOR = 2; // current
// Tail codec id.
inline constexpr uint8_t KPX_TAIL_VBYTE     = 1; // LEB128 varint stream (FOR base + values)

#pragma pack(push, 1)
struct KpxHeader {
    char     magic[4];                   // 0x00: "KPX4"
    uint16_t format_version;             // 0x04: 4
    uint8_t  k;                          // 0x06
    uint8_t  t;                          // 0x07: template length (0=contiguous)
    uint64_t total_postings;             // 0x08
    uint8_t  template_type;              // 0x10: TemplateType enum
    uint8_t  offset_type;                // 0x11: 0=uint32 offsets, 1=uint64 offsets
    uint8_t  reserved2[14];              // 0x12
    // ---- v4 codec extension (32 B, total header size = 64 B) ----
    uint8_t  codec_id;                   // 0x20: KPX_CODEC_PFOR_S8B
    uint8_t  codec_version;              // 0x21: 1
    uint16_t block_size;                 // 0x22: 128
    uint8_t  tail_codec;                 // 0x24: KPX_TAIL_VBYTE
    uint8_t  reserved0[27];              // 0x25
};
#pragma pack(pop)

static_assert(sizeof(KpxHeader) == 64, "KpxHeader v4 must be 64 bytes");

} // namespace ikafssn
