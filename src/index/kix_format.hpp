#pragma once

#include <cstdint>
#include <cstring>

namespace ikafssn {

// .kix v7 (Phase 5i): per-posting payload is FastPFor's
// CompositeCodec<SIMDFastPFor<4>, VariableByte> (PForDelta + VByte
// exception stream) — same codec as v6 — but the input stream is now the
// **distinct seq_id delta list** for the k-mer, i.e.
//   [abs_first, d1, d2, ...]   with d_i >= 1
// (intra-sequence k-mer duplicates are removed by a SIMD dedup kernel at
// build time).  As a result the leading u32 of every posting blob is the
// distinct seq_id count (renamed from `count` to `distinct_count`), and
// the header carries `total_distinct_postings` instead of `total_postings`.
inline constexpr char KIX_MAGIC[4] = {'K', 'I', 'X', '7'};

// Flag bits
inline constexpr uint32_t KIX_FLAG_SEQ_ID_WIDTH = 0x01; // 0=uint32, 1=uint64 (future)
inline constexpr uint32_t KIX_FLAG_HAS_KSX      = 0x02; // 0=no .ksx, 1=has .ksx
inline constexpr uint32_t KIX_FLAG_OFFSET32     = 0x04; // 0=uint64 offsets, 1=uint32 offsets

// Codec selection follows format_version since Phase 5g-1.  The legacy
// codec_id field in the header is kept as a reserved byte for header-size
// stability but readers do not dispatch on it.

#pragma pack(push, 1)
struct KixHeader {
    char     magic[4];                   // 0x00: "KIX7"
    uint16_t format_version;             // 0x04: 7
    uint8_t  k;                          // 0x06
    uint8_t  kmer_type;                  // 0x07: 0=uint16, 1=uint32
    uint32_t num_sequences;              // 0x08
    uint64_t total_distinct_postings;    // 0x0C: sum of distinct seq_id counts across all k-mers
    uint32_t flags;                      // 0x14
    uint16_t volume_index;               // 0x18
    uint16_t total_volumes;              // 0x1A
    uint16_t db_len;                     // 0x1C
    uint8_t  t;                          // 0x1E: template length (0=contiguous)
    uint8_t  template_type;              // 0x1F: TemplateType enum (0=contiguous)
    char     db[32];                     // 0x20
    // ---- 32 B codec-extension area (total header size = 96 B) ----
    // codec_id / codec_version / block_size / tail_codec are reserved
    // since Phase 5g-1: codec selection now follows format_version. They
    // are kept as struct fields for header-byte stability and so that
    // readers built against an older header give a clean format_version
    // mismatch error rather than misparsing the codec area.
    uint8_t  codec_id;                   // 0x40: 0 (reserved)
    uint8_t  codec_version;              // 0x41: 0 (reserved)
    uint16_t block_size;                 // 0x42: 0 (reserved)
    uint8_t  tail_codec;                 // 0x44: 0 (reserved)
    uint8_t  reserved0[3];               // 0x45
    uint32_t exception_codec_flags;      // 0x48: future use
    uint8_t  reserved1[20];              // 0x4C: pads header to a 96 B total
};
#pragma pack(pop)

static_assert(sizeof(KixHeader) == 96, "KixHeader v7 must be 96 bytes");

} // namespace ikafssn
