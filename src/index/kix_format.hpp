#pragma once

#include <cstdint>
#include <cstring>

namespace ikafssn {

// .kix v10 (Phase 1 of the fragment-indexing initiative): the wire layout
// of the dictionary and posting file bodies is unchanged from v9 (Elias-
// Fano dictionary + the v9 posting list shape — see kix_format v9 notes
// preserved in git history for the wire details).  v10 adds three header
// fields that drive short-sequence filtering and (in Phase 2 onward) long-
// sequence fragment splitting:
//
//   uint32_t min_seq_length    // -min_seq_length passed to ikafssnindex
//   uint32_t min_length_split  // 0 in Phase 1; non-zero from Phase 2
//   uint32_t overlap_length    // 0 in Phase 1; non-zero from Phase 2
//
// The magic is widened from 4 bytes "KIX9" to 5 bytes "KIX10"; the rest
// of the header shifts +1 and one byte is taken out of the trailing
// reserved tail so the total header size stays at 96 bytes.
inline constexpr char KIX_MAGIC[5] = {'K', 'I', 'X', '1', '0'};

// Flag bits
inline constexpr uint32_t KIX_FLAG_SEQ_ID_WIDTH = 0x01; // 0=uint32, 1=uint64 (future)
inline constexpr uint32_t KIX_FLAG_HAS_KSX      = 0x02; // 0=no .ksx, 1=has .ksx
inline constexpr uint32_t KIX_FLAG_OFFSET32     = 0x04; // 0=uint64 offsets, 1=uint32 offsets

// Codec selection follows format_version since Phase 5g-1.  The legacy
// codec_id field in the header is kept as a reserved byte for header-size
// stability but readers do not dispatch on it.

#pragma pack(push, 1)
struct KixHeader {
    char     magic[5];                   // 0x00: "KIX10"
    uint16_t format_version;             // 0x05: 10
    uint8_t  k;                          // 0x07
    uint8_t  kmer_type;                  // 0x08: 0=uint16, 1=uint32
    uint32_t num_sequences;              // 0x09: number of internal sequence ids (= fragments)
    uint64_t total_distinct_postings;    // 0x0D: sum of distinct seq_id counts across all k-mers
    uint32_t flags;                      // 0x15
    uint16_t volume_index;               // 0x19
    uint16_t total_volumes;              // 0x1B
    uint16_t db_len;                     // 0x1D
    uint8_t  t;                          // 0x1F: template length (0=contiguous)
    uint8_t  template_type;              // 0x20: TemplateType enum (0=contiguous)
    char     db[32];                     // 0x21
    // ---- 32 B codec-extension area (total header size = 96 B) ----
    // codec_id / codec_version / block_size / tail_codec are reserved
    // since Phase 5g-1: codec selection now follows format_version. They
    // are kept as struct fields for header-byte stability.
    uint8_t  codec_id;                   // 0x41: 0 (reserved)
    uint8_t  codec_version;              // 0x42: 0 (reserved)
    uint16_t block_size;                 // 0x43: 0 (reserved)
    uint8_t  tail_codec;                 // 0x45: 0 (reserved)
    uint8_t  reserved0[3];               // 0x46
    uint32_t exception_codec_flags;      // 0x49: future use
    // v10 additions (offset 0x4D, 12 bytes total)
    uint32_t min_seq_length;             // 0x4D: short-sequence filter threshold
    uint32_t min_length_split;           // 0x51: 0 in Phase 1
    uint32_t overlap_length;             // 0x55: 0 in Phase 1
    uint8_t  reserved1[7];               // 0x59: pads header to a 96 B total
};
#pragma pack(pop)

static_assert(sizeof(KixHeader) == 96, "KixHeader v10 must be 96 bytes");

} // namespace ikafssn
