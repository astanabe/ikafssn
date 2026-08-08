#include "core/kmer_revcomp_simd.hpp"
#include "core/kmer_encoding.hpp"
#include "util/simd_dispatch.hpp"

#include <cstddef>
#include <cstdint>
#include <type_traits>

#if defined(__x86_64__) || defined(__i386__)
  #include <immintrin.h>
  #define IKAFSSN_KR_X86 1
#endif

#if defined(__aarch64__)
  #include <arm_neon.h>
  #define IKAFSSN_KR_ARM 1
  #if __has_include(<arm_sve.h>)
    #include <arm_sve.h>
    #define IKAFSSN_KR_HAVE_SVE_HEADER 1
  #endif
#endif

namespace ikafssn {

// ===========================================================================
// Scalar reference: bit-exact source of truth.
// ===========================================================================
template <typename KmerInt>
static inline void kmer_revcomp_batch_scalar(const KmerInt* in, KmerInt* out,
                                             std::size_t n, int k) noexcept {
    for (std::size_t i = 0; i < n; ++i) out[i] = kmer_revcomp<KmerInt>(in[i], k);
}

// ===========================================================================
// x86_64 — uint32_t kernels (k = 9..16). Process 4 / 8 / 16 lanes per tier.
// All kernels share the same shape: ~v -> byte-reverse-within-dword (vpshufb)
// -> nibble swap -> 2-bit swap -> right shift to drop unused high bits.
// ===========================================================================
#if IKAFSSN_KR_X86

__attribute__((target("sse4.2,ssse3")))
static void kmer_revcomp_batch_sse42_u32(const std::uint32_t* in,
                                         std::uint32_t* out,
                                         std::size_t n, int k) noexcept {
    const int shift_out = 32 - 2 * k;
    const __m128i bswap = _mm_setr_epi8(
        3,2,1,0,  7,6,5,4,  11,10,9,8,  15,14,13,12);
    const __m128i mask_0f = _mm_set1_epi8(0x0F);
    const __m128i mask_33 = _mm_set1_epi8(0x33);
    std::size_t i = 0;
    for (; i + 4 <= n; i += 4) {
        __m128i v = _mm_loadu_si128(reinterpret_cast<const __m128i*>(in + i));
        v = _mm_xor_si128(v, _mm_set1_epi8(-1));
        v = _mm_shuffle_epi8(v, bswap);
        // Nibble swap within each byte: (v & 0x0F) << 4 | (v >> 4) & 0x0F
        __m128i lo = _mm_and_si128(v, mask_0f);
        __m128i hi = _mm_and_si128(_mm_srli_epi16(v, 4), mask_0f);
        v = _mm_or_si128(_mm_slli_epi16(lo, 4), hi);
        // 2-bit pair swap within each nibble
        __m128i lo2 = _mm_and_si128(v, mask_33);
        __m128i hi2 = _mm_and_si128(_mm_srli_epi16(v, 2), mask_33);
        v = _mm_or_si128(_mm_slli_epi16(lo2, 2), hi2);
        v = _mm_srli_epi32(v, shift_out);
        _mm_storeu_si128(reinterpret_cast<__m128i*>(out + i), v);
    }
    for (; i < n; ++i) out[i] = kmer_revcomp<std::uint32_t>(in[i], k);
}

__attribute__((target("avx2")))
static void kmer_revcomp_batch_avx2_u32(const std::uint32_t* in,
                                        std::uint32_t* out,
                                        std::size_t n, int k) noexcept {
    const int shift_out = 32 - 2 * k;
    const __m256i bswap = _mm256_setr_epi8(
        3,2,1,0,  7,6,5,4,  11,10,9,8,  15,14,13,12,
        3,2,1,0,  7,6,5,4,  11,10,9,8,  15,14,13,12);
    const __m256i mask_0f = _mm256_set1_epi8(0x0F);
    const __m256i mask_33 = _mm256_set1_epi8(0x33);
    std::size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        __m256i v = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(in + i));
        v = _mm256_xor_si256(v, _mm256_set1_epi8(-1));
        v = _mm256_shuffle_epi8(v, bswap);
        __m256i lo = _mm256_and_si256(v, mask_0f);
        __m256i hi = _mm256_and_si256(_mm256_srli_epi16(v, 4), mask_0f);
        v = _mm256_or_si256(_mm256_slli_epi16(lo, 4), hi);
        __m256i lo2 = _mm256_and_si256(v, mask_33);
        __m256i hi2 = _mm256_and_si256(_mm256_srli_epi16(v, 2), mask_33);
        v = _mm256_or_si256(_mm256_slli_epi16(lo2, 2), hi2);
        v = _mm256_srli_epi32(v, shift_out);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(out + i), v);
    }
    for (; i < n; ++i) out[i] = kmer_revcomp<std::uint32_t>(in[i], k);
}

__attribute__((target("avx512bw,avx512f,avx2")))
static void kmer_revcomp_batch_avx512bw_u32(const std::uint32_t* in,
                                            std::uint32_t* out,
                                            std::size_t n, int k) noexcept {
    const int shift_out = 32 - 2 * k;
    alignas(64) static const std::int8_t bswap_pattern[64] = {
        3,2,1,0,  7,6,5,4,  11,10,9,8,  15,14,13,12,
        3,2,1,0,  7,6,5,4,  11,10,9,8,  15,14,13,12,
        3,2,1,0,  7,6,5,4,  11,10,9,8,  15,14,13,12,
        3,2,1,0,  7,6,5,4,  11,10,9,8,  15,14,13,12,
    };
    const __m512i bswap = _mm512_load_si512(bswap_pattern);
    const __m512i mask_0f = _mm512_set1_epi8(0x0F);
    const __m512i mask_33 = _mm512_set1_epi8(0x33);
    std::size_t i = 0;
    for (; i + 16 <= n; i += 16) {
        __m512i v = _mm512_loadu_si512(reinterpret_cast<const void*>(in + i));
        v = _mm512_xor_si512(v, _mm512_set1_epi8(-1));
        v = _mm512_shuffle_epi8(v, bswap);
        __m512i lo = _mm512_and_si512(v, mask_0f);
        __m512i hi = _mm512_and_si512(_mm512_srli_epi16(v, 4), mask_0f);
        v = _mm512_or_si512(_mm512_slli_epi16(lo, 4), hi);
        __m512i lo2 = _mm512_and_si512(v, mask_33);
        __m512i hi2 = _mm512_and_si512(_mm512_srli_epi16(v, 2), mask_33);
        v = _mm512_or_si512(_mm512_slli_epi16(lo2, 2), hi2);
        v = _mm512_srli_epi32(v, shift_out);
        _mm512_storeu_si512(reinterpret_cast<void*>(out + i), v);
    }
    for (; i < n; ++i) out[i] = kmer_revcomp<std::uint32_t>(in[i], k);
}

__attribute__((target("avx512vbmi2,avx512vbmi,avx512bw,avx512f,avx2")))
static void kmer_revcomp_batch_avx512vbmi2_u32(const std::uint32_t* in,
                                               std::uint32_t* out,
                                               std::size_t n, int k) noexcept {
    // VBMI2 adds nothing this kernel can use; the body matches the VBMI/BW
    // tier.  Its own symbol so every SimdCap level dispatches explicitly.
    const int shift_out = 32 - 2 * k;
    alignas(64) static const std::int8_t bswap_pattern[64] = {
        3,2,1,0,  7,6,5,4,  11,10,9,8,  15,14,13,12,
        3,2,1,0,  7,6,5,4,  11,10,9,8,  15,14,13,12,
        3,2,1,0,  7,6,5,4,  11,10,9,8,  15,14,13,12,
        3,2,1,0,  7,6,5,4,  11,10,9,8,  15,14,13,12,
    };
    const __m512i bswap = _mm512_load_si512(bswap_pattern);
    const __m512i mask_0f = _mm512_set1_epi8(0x0F);
    const __m512i mask_33 = _mm512_set1_epi8(0x33);
    std::size_t i = 0;
    for (; i + 16 <= n; i += 16) {
        __m512i v = _mm512_loadu_si512(reinterpret_cast<const void*>(in + i));
        v = _mm512_xor_si512(v, _mm512_set1_epi8(-1));
        v = _mm512_shuffle_epi8(v, bswap);
        __m512i lo = _mm512_and_si512(v, mask_0f);
        __m512i hi = _mm512_and_si512(_mm512_srli_epi16(v, 4), mask_0f);
        v = _mm512_or_si512(_mm512_slli_epi16(lo, 4), hi);
        __m512i lo2 = _mm512_and_si512(v, mask_33);
        __m512i hi2 = _mm512_and_si512(_mm512_srli_epi16(v, 2), mask_33);
        v = _mm512_or_si512(_mm512_slli_epi16(lo2, 2), hi2);
        v = _mm512_srli_epi32(v, shift_out);
        _mm512_storeu_si512(reinterpret_cast<void*>(out + i), v);
    }
    for (; i < n; ++i) out[i] = kmer_revcomp<std::uint32_t>(in[i], k);
}

// ===========================================================================
// x86_64 — uint16_t kernels (k = 4..8). Same shape but per-u16 byte swap.
// ===========================================================================
__attribute__((target("sse4.2,ssse3")))
static void kmer_revcomp_batch_sse42_u16(const std::uint16_t* in,
                                         std::uint16_t* out,
                                         std::size_t n, int k) noexcept {
    const int shift_out = 16 - 2 * k;
    const __m128i bswap = _mm_setr_epi8(
        1,0, 3,2, 5,4, 7,6, 9,8, 11,10, 13,12, 15,14);
    const __m128i mask_0f = _mm_set1_epi8(0x0F);
    const __m128i mask_33 = _mm_set1_epi8(0x33);
    std::size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        __m128i v = _mm_loadu_si128(reinterpret_cast<const __m128i*>(in + i));
        v = _mm_xor_si128(v, _mm_set1_epi8(-1));
        v = _mm_shuffle_epi8(v, bswap);
        __m128i lo = _mm_and_si128(v, mask_0f);
        __m128i hi = _mm_and_si128(_mm_srli_epi16(v, 4), mask_0f);
        v = _mm_or_si128(_mm_slli_epi16(lo, 4), hi);
        __m128i lo2 = _mm_and_si128(v, mask_33);
        __m128i hi2 = _mm_and_si128(_mm_srli_epi16(v, 2), mask_33);
        v = _mm_or_si128(_mm_slli_epi16(lo2, 2), hi2);
        v = _mm_srli_epi16(v, shift_out);
        _mm_storeu_si128(reinterpret_cast<__m128i*>(out + i), v);
    }
    for (; i < n; ++i) out[i] = kmer_revcomp<std::uint16_t>(in[i], k);
}

__attribute__((target("avx2")))
static void kmer_revcomp_batch_avx2_u16(const std::uint16_t* in,
                                        std::uint16_t* out,
                                        std::size_t n, int k) noexcept {
    const int shift_out = 16 - 2 * k;
    const __m256i bswap = _mm256_setr_epi8(
        1,0, 3,2, 5,4, 7,6, 9,8, 11,10, 13,12, 15,14,
        1,0, 3,2, 5,4, 7,6, 9,8, 11,10, 13,12, 15,14);
    const __m256i mask_0f = _mm256_set1_epi8(0x0F);
    const __m256i mask_33 = _mm256_set1_epi8(0x33);
    std::size_t i = 0;
    for (; i + 16 <= n; i += 16) {
        __m256i v = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(in + i));
        v = _mm256_xor_si256(v, _mm256_set1_epi8(-1));
        v = _mm256_shuffle_epi8(v, bswap);
        __m256i lo = _mm256_and_si256(v, mask_0f);
        __m256i hi = _mm256_and_si256(_mm256_srli_epi16(v, 4), mask_0f);
        v = _mm256_or_si256(_mm256_slli_epi16(lo, 4), hi);
        __m256i lo2 = _mm256_and_si256(v, mask_33);
        __m256i hi2 = _mm256_and_si256(_mm256_srli_epi16(v, 2), mask_33);
        v = _mm256_or_si256(_mm256_slli_epi16(lo2, 2), hi2);
        v = _mm256_srli_epi16(v, shift_out);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(out + i), v);
    }
    for (; i < n; ++i) out[i] = kmer_revcomp<std::uint16_t>(in[i], k);
}

__attribute__((target("avx512bw,avx512f,avx2")))
static void kmer_revcomp_batch_avx512bw_u16(const std::uint16_t* in,
                                            std::uint16_t* out,
                                            std::size_t n, int k) noexcept {
    const int shift_out = 16 - 2 * k;
    alignas(64) static const std::int8_t bswap_pattern[64] = {
        1,0, 3,2, 5,4, 7,6, 9,8, 11,10, 13,12, 15,14,
        1,0, 3,2, 5,4, 7,6, 9,8, 11,10, 13,12, 15,14,
        1,0, 3,2, 5,4, 7,6, 9,8, 11,10, 13,12, 15,14,
        1,0, 3,2, 5,4, 7,6, 9,8, 11,10, 13,12, 15,14,
    };
    const __m512i bswap = _mm512_load_si512(bswap_pattern);
    const __m512i mask_0f = _mm512_set1_epi8(0x0F);
    const __m512i mask_33 = _mm512_set1_epi8(0x33);
    std::size_t i = 0;
    for (; i + 32 <= n; i += 32) {
        __m512i v = _mm512_loadu_si512(reinterpret_cast<const void*>(in + i));
        v = _mm512_xor_si512(v, _mm512_set1_epi8(-1));
        v = _mm512_shuffle_epi8(v, bswap);
        __m512i lo = _mm512_and_si512(v, mask_0f);
        __m512i hi = _mm512_and_si512(_mm512_srli_epi16(v, 4), mask_0f);
        v = _mm512_or_si512(_mm512_slli_epi16(lo, 4), hi);
        __m512i lo2 = _mm512_and_si512(v, mask_33);
        __m512i hi2 = _mm512_and_si512(_mm512_srli_epi16(v, 2), mask_33);
        v = _mm512_or_si512(_mm512_slli_epi16(lo2, 2), hi2);
        v = _mm512_srli_epi16(v, shift_out);
        _mm512_storeu_si512(reinterpret_cast<void*>(out + i), v);
    }
    for (; i < n; ++i) out[i] = kmer_revcomp<std::uint16_t>(in[i], k);
}

__attribute__((target("avx512vbmi2,avx512vbmi,avx512bw,avx512f,avx2")))
static void kmer_revcomp_batch_avx512vbmi2_u16(const std::uint16_t* in,
                                               std::uint16_t* out,
                                               std::size_t n, int k) noexcept {
    // Body matches BW; separate symbol for tier-effect measurement.
    kmer_revcomp_batch_avx512bw_u16(in, out, n, k);
}

#endif // IKAFSSN_KR_X86

// ===========================================================================
// aarch64 kernels.
// ===========================================================================
#if IKAFSSN_KR_ARM

static void kmer_revcomp_batch_neon_u32(const std::uint32_t* in,
                                        std::uint32_t* out,
                                        std::size_t n, int k) noexcept {
    const int shift_out = 32 - 2 * k;
    std::size_t i = 0;
    const uint8x16_t mask_0f = vdupq_n_u8(0x0F);
    const uint8x16_t mask_33 = vdupq_n_u8(0x33);
    for (; i + 4 <= n; i += 4) {
        uint8x16_t v = vld1q_u8(reinterpret_cast<const std::uint8_t*>(in + i));
        v = vmvnq_u8(v);
        // Byte-reverse within each 32-bit lane.
        v = vrev32q_u8(v);
        // Nibble swap
        uint8x16_t lo = vandq_u8(v, mask_0f);
        uint8x16_t hi = vandq_u8(vshrq_n_u8(v, 4), mask_0f);
        v = vorrq_u8(vshlq_n_u8(lo, 4), hi);
        // 2-bit pair swap
        uint8x16_t lo2 = vandq_u8(v, mask_33);
        uint8x16_t hi2 = vandq_u8(vshrq_n_u8(v, 2), mask_33);
        v = vorrq_u8(vshlq_n_u8(lo2, 2), hi2);
        // Right shift each u32 by shift_out
        uint32x4_t v32 = vreinterpretq_u32_u8(v);
        v32 = vshlq_u32(v32, vdupq_n_s32(-shift_out));
        vst1q_u32(out + i, v32);
    }
    for (; i < n; ++i) out[i] = kmer_revcomp<std::uint32_t>(in[i], k);
}

static void kmer_revcomp_batch_neon_u16(const std::uint16_t* in,
                                        std::uint16_t* out,
                                        std::size_t n, int k) noexcept {
    const int shift_out = 16 - 2 * k;
    std::size_t i = 0;
    const uint8x16_t mask_0f = vdupq_n_u8(0x0F);
    const uint8x16_t mask_33 = vdupq_n_u8(0x33);
    for (; i + 8 <= n; i += 8) {
        uint8x16_t v = vld1q_u8(reinterpret_cast<const std::uint8_t*>(in + i));
        v = vmvnq_u8(v);
        v = vrev16q_u8(v);
        uint8x16_t lo = vandq_u8(v, mask_0f);
        uint8x16_t hi = vandq_u8(vshrq_n_u8(v, 4), mask_0f);
        v = vorrq_u8(vshlq_n_u8(lo, 4), hi);
        uint8x16_t lo2 = vandq_u8(v, mask_33);
        uint8x16_t hi2 = vandq_u8(vshrq_n_u8(v, 2), mask_33);
        v = vorrq_u8(vshlq_n_u8(lo2, 2), hi2);
        uint16x8_t v16 = vreinterpretq_u16_u8(v);
        v16 = vshlq_u16(v16, vdupq_n_s16(static_cast<int16_t>(-shift_out)));
        vst1q_u16(out + i, v16);
    }
    for (; i < n; ++i) out[i] = kmer_revcomp<std::uint16_t>(in[i], k);
}

#if IKAFSSN_KR_HAVE_SVE_HEADER

__attribute__((target("+sve")))
static void kmer_revcomp_batch_sve_u32(const std::uint32_t* in,
                                       std::uint32_t* out,
                                       std::size_t n, int k) noexcept {
    // SVE has no single-instruction byte-reverse-per-dword (vrev32q_u8), and
    // SVE/SVE2 gains on this pattern are small, so a scalar body inside an
    // SVE-attribute function is the simplest correct path.
    (void)k;
    for (std::size_t i = 0; i < n; ++i) out[i] = kmer_revcomp<std::uint32_t>(in[i], k);
}

__attribute__((target("+sve2")))
static void kmer_revcomp_batch_sve2_u32(const std::uint32_t* in,
                                        std::uint32_t* out,
                                        std::size_t n, int k) noexcept {
    for (std::size_t i = 0; i < n; ++i) out[i] = kmer_revcomp<std::uint32_t>(in[i], k);
}

__attribute__((target("+sve")))
static void kmer_revcomp_batch_sve_u16(const std::uint16_t* in,
                                       std::uint16_t* out,
                                       std::size_t n, int k) noexcept {
    for (std::size_t i = 0; i < n; ++i) out[i] = kmer_revcomp<std::uint16_t>(in[i], k);
}

__attribute__((target("+sve2")))
static void kmer_revcomp_batch_sve2_u16(const std::uint16_t* in,
                                        std::uint16_t* out,
                                        std::size_t n, int k) noexcept {
    for (std::size_t i = 0; i < n; ++i) out[i] = kmer_revcomp<std::uint16_t>(in[i], k);
}

#endif // IKAFSSN_KR_HAVE_SVE_HEADER

#endif // IKAFSSN_KR_ARM

// ===========================================================================
// Public dispatcher.
// ===========================================================================
template <typename KmerInt>
void kmer_revcomp_batch(const KmerInt* in, KmerInt* out, std::size_t n,
                        int k) noexcept {
    if (n == 0) return;
    SimdCap cap = current_simd_cap();

#if IKAFSSN_KR_X86
    if constexpr (std::is_same_v<KmerInt, std::uint32_t>) {
        if (cap >= SimdCap::AVX512VBMI2 && n >= 16) {
            kmer_revcomp_batch_avx512vbmi2_u32(in, out, n, k); return;
        }
        if (cap >= SimdCap::AVX512BW && n >= 16) {
            kmer_revcomp_batch_avx512bw_u32(in, out, n, k); return;
        }
        if (cap >= SimdCap::AVX2 && n >= 8) {
            kmer_revcomp_batch_avx2_u32(in, out, n, k); return;
        }
        if (cap >= SimdCap::SSE42 && n >= 4) {
            kmer_revcomp_batch_sse42_u32(in, out, n, k); return;
        }
    } else if constexpr (std::is_same_v<KmerInt, std::uint16_t>) {
        if (cap >= SimdCap::AVX512VBMI2 && n >= 32) {
            kmer_revcomp_batch_avx512vbmi2_u16(in, out, n, k); return;
        }
        if (cap >= SimdCap::AVX512BW && n >= 32) {
            kmer_revcomp_batch_avx512bw_u16(in, out, n, k); return;
        }
        if (cap >= SimdCap::AVX2 && n >= 16) {
            kmer_revcomp_batch_avx2_u16(in, out, n, k); return;
        }
        if (cap >= SimdCap::SSE42 && n >= 8) {
            kmer_revcomp_batch_sse42_u16(in, out, n, k); return;
        }
    }
#elif IKAFSSN_KR_ARM
    if constexpr (std::is_same_v<KmerInt, std::uint32_t>) {
      #if IKAFSSN_KR_HAVE_SVE_HEADER
        if (cap >= SimdCap::SVE2 && n >= 8) {
            kmer_revcomp_batch_sve2_u32(in, out, n, k); return;
        }
        if (cap >= SimdCap::SVE && n >= 8) {
            kmer_revcomp_batch_sve_u32(in, out, n, k); return;
        }
      #endif
        if (cap >= SimdCap::NEON && n >= 4) {
            kmer_revcomp_batch_neon_u32(in, out, n, k); return;
        }
    } else if constexpr (std::is_same_v<KmerInt, std::uint16_t>) {
      #if IKAFSSN_KR_HAVE_SVE_HEADER
        if (cap >= SimdCap::SVE2 && n >= 16) {
            kmer_revcomp_batch_sve2_u16(in, out, n, k); return;
        }
        if (cap >= SimdCap::SVE && n >= 16) {
            kmer_revcomp_batch_sve_u16(in, out, n, k); return;
        }
      #endif
        if (cap >= SimdCap::NEON && n >= 8) {
            kmer_revcomp_batch_neon_u16(in, out, n, k); return;
        }
    }
#else
    (void)cap;
#endif
    kmer_revcomp_batch_scalar<KmerInt>(in, out, n, k);
}

template void kmer_revcomp_batch<std::uint16_t>(
    const std::uint16_t*, std::uint16_t*, std::size_t, int) noexcept;
template void kmer_revcomp_batch<std::uint32_t>(
    const std::uint32_t*, std::uint32_t*, std::size_t, int) noexcept;

} // namespace ikafssn
