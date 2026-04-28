#pragma once

#include <cstdint>

namespace ikafssn {

// .kpx v6 (Phase 5g-2): per-(kmer, seq_id) partitioned position posting.
// Each (k-mer, seq_id) cluster whose occurrence count exceeds the build-
// time `freq_threshold_part` is split out into its own partition group;
// the rest of the occurrences for that k-mer are merged into a single
// short bucket.  Both partition groups and the short bucket use the
// Phase 5e FOR-within-block stream layout.  See src/index/pfd_codec.hpp
// for the byte-level wire format and lock-step merge protocol.
inline constexpr char KPX_MAGIC[4] = {'K', 'P', 'X', '6'};

// Codec selection follows format_version (since Phase 5g-1). The header
// still carries an 8-bit codec_id field for layout stability but readers
// no longer dispatch on it.

#pragma pack(push, 1)
struct KpxHeader {
    char     magic[4];                   // 0x00: "KPX6"
    uint16_t format_version;             // 0x04: 6
    uint8_t  k;                          // 0x06
    uint8_t  t;                          // 0x07: template length (0=contiguous)
    uint64_t total_postings;             // 0x08
    uint8_t  template_type;              // 0x10: TemplateType enum
    uint8_t  offset_type;                // 0x11: 0=uint32 offsets, 1=uint64 offsets
    uint8_t  reserved2[14];              // 0x12
    // ---- 32 B codec-extension area (total header size = 64 B) ----
    // codec_id / codec_version / block_size / tail_codec are reserved
    // since Phase 5g-1: codec selection now follows format_version. Kept
    // as fields for header-byte stability.
    uint8_t  codec_id;                   // 0x20: 0 (reserved)
    uint8_t  codec_version;              // 0x21: 0 (reserved)
    uint16_t block_size;                 // 0x22: 0 (reserved)
    uint8_t  tail_codec;                 // 0x24: 0 (reserved)
    uint8_t  reserved0[27];              // 0x25
};
#pragma pack(pop)

static_assert(sizeof(KpxHeader) == 64, "KpxHeader v6 must be 64 bytes");

} // namespace ikafssn
