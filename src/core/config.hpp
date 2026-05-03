#pragma once

#include <cstdint>
#include <cstddef>

#if __cplusplus >= 202002L
#include <bit>
#endif

namespace ikafssn {

// Endianness check
#if __cplusplus >= 202002L
static_assert(std::endian::native == std::endian::little,
              "ikafssn requires a little-endian platform");
#else
static_assert(__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__,
              "ikafssn requires a little-endian platform");
#endif

// k-mer length limits
//
// MAX_K = 15 (Phase 7-post): the previous nominal upper of 16 was
// unreachable in practice because table_size(16) = 4^16 = 2^32
// overflows the uint32_t return type and silently wraps to 0.  k=15
// (table_size = 4^15 = 1,073,741,824) is the largest k that fits in
// every uint32_t array index this codebase uses (kix dictionary
// length, .khx bitset, EFHeader.D / upper_bits).  k > 15 is also
// unattractive in practice — the resulting dictionary blob and
// posting file blow up the per-volume RAM budget without commensurate
// search-quality gain (use minimizer-style indexes for that regime).
//
// KmerScanner / packed_kmer_scanner internally hold MAX_K-sized
// arrays (window_degen, AmbigInfo[]) so growing MAX_K back to 16
// would also require widening those buffers.
inline constexpr int MIN_K = 5;
inline constexpr int MAX_K = 15;
inline constexpr int K_TYPE_THRESHOLD = 9; // k >= 9 uses uint32_t (contiguous/t=0 only; for spaced seeds use kmer_type_for(k, t))

// Format versions — all index files share a single major version so users
// only have to track one number. Phase 5c established this convention at
// v4. Phase 5g-1 bumped to v5 (.kix codec replaced with FastPFor's
// CompositeCodec<SIMDFastPFor<4>, VariableByte>). Phase 5g-2 bumped to v6
// (.kpx partition + short bucket layout). Phase 5i bumped to v7 (distinct
// seq_id semantic + self-describing short bucket). Phase 6 bumped to v8
// (kind-map .kpx layout). Phase 7 bumps to v9:
//   - .kix and .kpx dictionaries replaced by Elias-Fano blobs (4.6× smaller
//     on average across NCBI nt_v4).  Per-tier ladder mirrors the FastPFor
//     PForDelta dispatcher: AVX2+ uses BMI2 PDEP+TZCNT for inner-word
//     select1; SSE4.2 / NEON use the popcount + bit-stripping fallback.
//     KIX_FLAG_OFFSET32 / KpxHeader.offset_type are reserved (writers force
//     the EF sentinel; readers ignore the byte).
//   - .kix posting list header dropped [u32 body_words] (Phase 7c dedup B);
//     body_words is derived from the EF dictionary's posting_byte_length.
//   - .kpx posting list header dropped all four redundant u32 fields
//     ([u32 distinct_count] (dedup A — comes from kix_count), the three
//     per-kind counts (dedup C — popcount of the 2-bit kind map), and
//     [u32 short2_position_count] (dedup D — horizontal sum of the u8
//     occ_count[] array)).  The posting list now starts directly at the
//     2-bit kind map; empty .kpx posting lists emit 0 bytes.
//   - .ksx / .khx data layouts unchanged; format_version bumps for
//     family-wide alignment.  .kvx FORMAT_VERSION line bumps to 9.
inline constexpr uint16_t KIX_FORMAT_VERSION = 9;
inline constexpr uint16_t KPX_FORMAT_VERSION = 9;
inline constexpr uint16_t KSX_FORMAT_VERSION = 9;
inline constexpr uint16_t KHX_FORMAT_VERSION = 9;

// Direct-address table size for k-mer value k: 4^k.
// Max supported: 4^15 = 1,073,741,824 (fits uint32_t).  k=16 would
// produce 4^16 = 2^32 which overflows the uint32_t return type and
// wraps to 0 — see MAX_K above.
inline constexpr uint32_t table_size(int k) {
    return static_cast<uint32_t>(uint64_t(1) << (2 * k));
}

// Mask for k-mer of given k: (1 << 2k) - 1
template <typename KmerInt>
inline constexpr KmerInt kmer_mask(int k) {
    return static_cast<KmerInt>((uint64_t(1) << (2 * k)) - 1);
}

// Align value up to given alignment (must be power of 2)
inline constexpr uint64_t align_up(uint64_t value, uint64_t alignment) {
    return (value + alignment - 1) & ~(alignment - 1);
}

} // namespace ikafssn
