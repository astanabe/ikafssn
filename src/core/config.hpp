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
inline constexpr int K_TYPE_THRESHOLD = 9; // k >= 9 uses uint32_t

// Format versions
inline constexpr uint16_t KIX_FORMAT_VERSION = 1;
inline constexpr uint16_t KPX_FORMAT_VERSION = 1;
inline constexpr uint16_t KSX_FORMAT_VERSION = 1;
inline constexpr uint16_t KHX_FORMAT_VERSION = 1;

// Direct-address table size for k-mer value k: 4^k
inline constexpr uint64_t table_size(int k) {
    return uint64_t(1) << (2 * k);
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
