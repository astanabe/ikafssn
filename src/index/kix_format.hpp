#pragma once

#include <cstdint>
#include <cstring>

namespace ikafssn {

// .kix v8 (Phase 6): on-wire payload layout is identical to v7 (per-posting
// FastPFor CompositeCodec<SIMDFastPFor<4>, VariableByte> over the
// distinct seq_id delta stream [abs_first, d1, d2, ...] with d_i >= 1).
// The format_version + magic bump aligns .kix with the v8 .kpx codec
// rewrite that occurs in the same phase (see src/index/kpx_format.hpp);
// the .kix decoded distinct_seq_id array is now part of the .kpx
// candidate-resolution contract and must always be available alongside
// the .kpx posting.
inline constexpr char KIX_MAGIC[4] = {'K', 'I', 'X', '8'};

// Flag bits
inline constexpr uint32_t KIX_FLAG_SEQ_ID_WIDTH = 0x01; // 0=uint32, 1=uint64 (future)
inline constexpr uint32_t KIX_FLAG_HAS_KSX      = 0x02; // 0=no .ksx, 1=has .ksx
inline constexpr uint32_t KIX_FLAG_OFFSET32     = 0x04; // 0=uint64 offsets, 1=uint32 offsets

// Codec selection follows format_version since Phase 5g-1.  The legacy
// codec_id field in the header is kept as a reserved byte for header-size
// stability but readers do not dispatch on it.

#pragma pack(push, 1)
struct KixHeader {
    char     magic[4];                   // 0x00: "KIX8"
    uint16_t format_version;             // 0x04: 8
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

static_assert(sizeof(KixHeader) == 96, "KixHeader v8 must be 96 bytes");

} // namespace ikafssn
