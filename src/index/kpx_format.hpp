#pragma once

#include <cstdint>

namespace ikafssn {

// .kpx v9 (Phase 7): per-(kmer, seq_id) partitioned position posting file.
// The pos_offsets dictionary is now an Elias-Fano blob (replacing the
// raw u32/u64 offsets — see src/index/ef_codec.hpp), and all four
// redundant per-posting-list header u32 fields were dropped (Phase 7c
// dedup A + 7d dedup C/D); the body starts directly at the 2-bit kind map.
// Each distinct seq_id is classified into one of three kinds via the
// kind map indexed by the .kix decoded distinct_seq_id array (the seq_id
// is therefore not stored in the .kpx posting list at all):
//
//   00 = short_occ1     — exactly 1 position
//   01 = short_occ_ge2  — between 2 and freq_threshold_part positions
//   10 = partition      — strictly more than freq_threshold_part positions
//   11 = reserved
//
// Per-posting-list layout (no fixed header — the body starts at the
// kind map):
//   [2-bit kind map: ceil(distinct_count*2/8) bytes]
//                                 // distinct_count comes from kix_count
//   for each partition (in .kix sid order):
//     [u32 occ_count][FOR stream over occ_count positions]
//   [FOR stream over short1_count positions]
//                                 // short1_count = popcount of kind map's
//                                 // 00 pairs
//   [u8  occ_count[short2_count]] // 2..freq_threshold_part
//                                 // short2_count = popcount of 01 pairs
//   [FOR stream over short2_position_count positions]
//                                 // short2_position_count = horizontal
//                                 // sum of u8 occ_count[] (or equivalently
//                                 // the running cum that builds the
//                                 // offset table)
//
// Empty .kpx posting lists (kix_count == 0) emit 0 bytes.
//
// Each FOR-block header is 8 B (proposal F): [u32 min][u8 b][3 B pad]
// followed by 16*b bytes bitpacked (value - min).  Stream tails switch
// from a varint stream to a packed bit-width stream (proposal D):
//   [u8 tail_count][u32 tail_min][u8 tail_b][bitpacked: ceil(tail_count*tail_b/8) B]
// (when tail_count == 0 only the leading u8 0 is written).
//
// Decoding is a 2-pointer merge walk of the .kix decoded distinct_seq_id
// array against the (sorted) candidate set: the kind map gives each
// candidate's (kind, rank) pair, which then indexes into the decoded
// partition / short1 / short2 sub-buckets.  See src/index/pfd_codec.hpp
// for the byte-level wire format and the open_stream_kpx_for_candidates
// API contract.
inline constexpr char KPX_MAGIC[4] = {'K', 'P', 'X', '9'};

// Codec selection follows format_version (since Phase 5g-1). The header
// still carries an 8-bit codec_id field for layout stability but readers
// no longer dispatch on it.

#pragma pack(push, 1)
struct KpxHeader {
    char     magic[4];                   // 0x00: "KPX9"
    uint16_t format_version;             // 0x04: 9
    uint8_t  k;                          // 0x06
    uint8_t  t;                          // 0x07: template length (0=contiguous)
    uint64_t total_position_count;       // 0x08: sum of position counts across all k-mers
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

static_assert(sizeof(KpxHeader) == 64, "KpxHeader v8 must be 64 bytes");

} // namespace ikafssn
