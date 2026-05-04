// Phase 7b — Elias-Fano dictionary codec, per-tier SIMD body.
//
// This translation unit is compiled once per ISA tier via the
// ikafssn_ef_<tier> OBJECT libraries declared in the top-level
// CMakeLists.txt.  Each compilation receives:
//
//   -DIKAFSSN_EF_TIER_NAME=<tier>   names the per-tier inner namespace
//                                    (ikafssn_ef_sse42 / ikafssn_ef_avx2 /
//                                     ikafssn_ef_avx512bw /
//                                     ikafssn_ef_avx512vbmi2 / ikafssn_ef_neon)
//   -m<arch> ...                     controls the instructions the
//                                    compiler may emit
//
// The hot path is the EF random-access ``access(i)``: locate the i-th
// set bit in the upper-bits unary stream via ``select1_after``, then
// concatenate the corresponding low-l bits.  When the build sees BMI2
// (AVX2 tier and above on x86; not present on SSE4.2 or NEON) the
// inner-word select uses ``PDEP`` + ``TZCNT`` for an O(1) per-word
// answer; otherwise we fall back to ``__builtin_popcountll`` +
// ``__builtin_ctzll`` with a tight bit-stripping loop.

#include "index/ef_codec.hpp"
#include "index/ef_format.hpp"

#include <cstdint>
#include <cstddef>
#include <cstring>
#include <vector>

#ifndef IKAFSSN_EF_TIER_NAME
#error "IKAFSSN_EF_TIER_NAME must be set per tier"
#endif

#if defined(__BMI2__)
  #include <immintrin.h>
#endif

#define IKAFSSN_EF_TIER_NS_(x) ikafssn_ef_##x
#define IKAFSSN_EF_TIER_NS(x)  IKAFSSN_EF_TIER_NS_(x)

namespace ikafssn::ef::IKAFSSN_EF_TIER_NS(IKAFSSN_EF_TIER_NAME) {

namespace {

inline std::uint64_t mask_for_bits(std::uint8_t l) noexcept {
    return (l >= 64) ? ~std::uint64_t{0}
                     : ((std::uint64_t{1} << l) - std::uint64_t{1});
}

inline std::uint64_t read_lower_bits(const std::uint64_t* lower,
                                     std::uint64_t bit_pos,
                                     std::uint8_t l,
                                     std::uint64_t mask_l) noexcept {
    if (l == 0) return 0;
    std::uint64_t word_idx    = bit_pos >> 6;
    std::uint64_t bit_in_word = bit_pos & 63;
    std::uint64_t lo = (lower[word_idx] >> bit_in_word) & mask_l;
    if (bit_in_word + l > 64) {
        lo |= (lower[word_idx + 1] << (64 - bit_in_word)) & mask_l;
    }
    return lo;
}

// Inner-word select: position (0..63) of the i-th set bit in `w`.  i is
// 0-indexed; the caller guarantees __builtin_popcountll(w) > i.  On
// AVX2-and-above tiers BMI2's PDEP+TZCNT collapses to two cycles; the
// SSE4.2 / NEON fallback strips i low bits then ctz.
inline std::uint64_t select_in_word(std::uint64_t w, std::uint32_t i) noexcept {
#if defined(__BMI2__)
    return static_cast<std::uint64_t>(
        _tzcnt_u64(_pdep_u64(std::uint64_t{1} << i, w)));
#else
    for (std::uint32_t k = 0; k < i; ++k) {
        w &= w - 1;
    }
    return static_cast<std::uint64_t>(__builtin_ctzll(w));
#endif
}

// Locate the bit position of the j-th (0-indexed) set bit AT OR AFTER
// upper-bit position `start_bit_pos`.  Word-by-word popcount + select
// in-word.
inline std::uint64_t select1_after(const std::uint64_t* upper,
                                   std::uint64_t upper_bits_total,
                                   std::uint64_t start_bit_pos,
                                   std::uint64_t j) noexcept {
    std::uint64_t pos = start_bit_pos;
    std::uint64_t found = 0;
    while (pos < upper_bits_total) {
        std::uint64_t word_idx    = pos >> 6;
        std::uint64_t bit_in_word = pos & 63;
        std::uint64_t w = upper[word_idx] >> bit_in_word;
        std::uint64_t cnt = static_cast<std::uint64_t>(__builtin_popcountll(w));
        if (found + cnt > j) {
            std::uint32_t need = static_cast<std::uint32_t>(j - found);
            return pos + select_in_word(w, need);
        }
        found += cnt;
        pos = (word_idx + 1) << 6;
    }
    return upper_bits_total;  // not found (caller error)
}

} // anonymous namespace

std::uint64_t access_dictionary_ef(const EFHeader& hdr,
                                   const std::uint64_t* lower,
                                   const std::uint64_t* upper,
                                   const std::uint64_t* select,
                                   std::uint32_t i) noexcept {
    const std::uint8_t l = hdr.l;
    const std::uint64_t mask_l = mask_for_bits(l);
    const std::uint64_t upper_bits_total = static_cast<std::uint64_t>(hdr.upper_bits);

    std::uint32_t sample_idx = i / kSelectStep;
    std::uint32_t base_rank  = sample_idx * kSelectStep;
    std::uint64_t base_pos   = select[sample_idx];
    std::uint64_t up_pos;
    if (i == base_rank) {
        up_pos = base_pos;
    } else {
        up_pos = select1_after(upper,
                               upper_bits_total,
                               base_pos + 1,
                               static_cast<std::uint64_t>(i - base_rank - 1));
    }

    std::uint64_t hi = up_pos - static_cast<std::uint64_t>(i);
    std::uint64_t lo = read_lower_bits(lower,
                                       static_cast<std::uint64_t>(i) * l,
                                       l, mask_l);
    std::uint64_t ef_value = (hi << l) | lo;
    return ef_value - static_cast<std::uint64_t>(i);
}

} // namespace ikafssn::ef::ikafssn_ef_<tier>
