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
inline constexpr int MIN_K = 5;
inline constexpr int MAX_K = 16;
inline constexpr int K_TYPE_THRESHOLD = 9; // k >= 9 uses uint32_t (contiguous/t=0 only; for spaced seeds use kmer_type_for(k, t))

// Format versions — all index files share a single major version so users
// only have to track one number. Phase 5c established this convention at
// v4. Phase 5g-1 bumped to v5 (.kix codec replaced with FastPFor's
// CompositeCodec<SIMDFastPFor<4>, VariableByte>). Phase 5g-2 bumped to v6
// (.kpx partition + short bucket layout). Phase 5i bumped to v7 (distinct
// seq_id semantic + self-describing short bucket). Phase 6 bumps to v8:
//   - .kpx short bucket is split into occ=1 (positions only, no
//     occ_count) and occ>=2 (u8 occ_count + positions) sub-buckets.
//   - .kpx classifies each distinct seq_id via a 2-bit kind map indexed
//     by the .kix decoded distinct_seq_id array; the per-cluster seq_id
//     duplicate is dropped from both partition groups and the short
//     buckets (only the .kix posting carries seq_ids).
//   - FOR-block header layout becomes [u32 min][u8 b][3 B pad] (8 B,
//     body 8 B aligned).
//   - FOR-block tail switches from a varint stream to a packed bit-width
//     stream: [u8 tail_count][u32 tail_min][u8 tail_b][bitpacked body].
//   - .kix / .ksx / .khx / .kvx data layouts unchanged; format_version
//     bumps for family-wide alignment (Phase 5c policy).
inline constexpr uint16_t KIX_FORMAT_VERSION = 8;
inline constexpr uint16_t KPX_FORMAT_VERSION = 8;
inline constexpr uint16_t KSX_FORMAT_VERSION = 8;
inline constexpr uint16_t KHX_FORMAT_VERSION = 8;

// Direct-address table size for k-mer value k: 4^k
// Max supported: 2 * 4^12 = 33,554,432 (fits uint32_t).
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
