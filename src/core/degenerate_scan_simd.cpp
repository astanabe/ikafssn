#include "core/degenerate_scan_simd.hpp"
#include "core/kmer_encoding.hpp"
#include "util/simd_dispatch.hpp"

#include <cstddef>
#include <cstdint>

#if defined(__x86_64__) || defined(__i386__)
  #include <immintrin.h>
  #define IKAFSSN_DS_X86 1
#endif

#if defined(__aarch64__)
  #include <arm_neon.h>
  #define IKAFSSN_DS_ARM 1
#endif

namespace ikafssn {

// Two-LUT set membership (Mula style).
//
// Degenerate IUPAC chars and their hex codes:
//   uppercase: B(0x42) D(0x44) H(0x48) K(0x4B) M(0x4D) N(0x4E)
//              R(0x52) S(0x53) V(0x56) W(0x57) Y(0x59)
//   lowercase: b(0x62) d(0x64) h(0x68) k(0x6B) m(0x6D) n(0x6E)
//              r(0x72) s(0x73) v(0x76) w(0x77) y(0x79)
//
// LUT_LO_BITSET[lo nibble] sets one bit per row that contains that low nibble:
//   bit 4 = uppercase row 0x4 (B,D,H,K,M,N)
//   bit 5 = uppercase row 0x5 (R,S,V,W,Y)
//   bit 6 = lowercase row 0x6 (b,d,h,k,m,n)
//   bit 7 = lowercase row 0x7 (r,s,v,w,y)
// LUT_HI_ROWMASK[hi nibble] selects which row applies (0x10/0x20/0x40/0x80
// for hi=4/5/6/7, else 0). Per-byte AND of the two LUT lookups is non-zero
// iff the byte is one of the 22 degenerate codes.

alignas(16) static const std::int8_t kLutLoBitset[16] = {
    0x00,        // 0x0
    0x00,        // 0x1
    static_cast<std::int8_t>(0xF0),  // 0x2: B,R,b,r -> bits 4,5,6,7
    static_cast<std::int8_t>(0xA0),  // 0x3: S,s     -> bits 5,7
    0x50,        // 0x4: D,d     -> bits 4,6
    0x00,        // 0x5
    static_cast<std::int8_t>(0xA0),  // 0x6: V,v     -> bits 5,7
    static_cast<std::int8_t>(0xA0),  // 0x7: W,w     -> bits 5,7
    0x50,        // 0x8: H,h     -> bits 4,6
    static_cast<std::int8_t>(0xA0),  // 0x9: Y,y     -> bits 5,7
    0x00,        // 0xA
    0x50,        // 0xB: K,k     -> bits 4,6
    0x00,        // 0xC
    0x50,        // 0xD: M,m     -> bits 4,6
    0x50,        // 0xE: N,n     -> bits 4,6
    0x00,        // 0xF
};

alignas(16) static const std::int8_t kLutHiRowmask[16] = {
    0x00, 0x00, 0x00, 0x00,
    0x10, 0x20, 0x40, static_cast<std::int8_t>(0x80),  // hi=4..7
    0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00,
};

// Scalar reference (always compiled): tail and SIMD-disabled fallback.
static inline bool has_degenerate_base_scalar(const char* data,
                                              std::size_t n) noexcept {
    const bool* tbl = degenerate_base_table();
    for (std::size_t i = 0; i < n; ++i) {
        if (tbl[static_cast<std::uint8_t>(data[i])]) return true;
    }
    return false;
}

// x86_64 kernels.
#if IKAFSSN_DS_X86

__attribute__((target("sse4.2,ssse3")))
static bool has_degenerate_base_sse42(const char* data,
                                      std::size_t n) noexcept {
    const __m128i lut_lo = _mm_load_si128(
        reinterpret_cast<const __m128i*>(kLutLoBitset));
    const __m128i lut_hi = _mm_load_si128(
        reinterpret_cast<const __m128i*>(kLutHiRowmask));
    const __m128i mask_0f = _mm_set1_epi8(0x0F);
    std::size_t i = 0;
    for (; i + 16 <= n; i += 16) {
        __m128i v  = _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(data + i));
        __m128i lo = _mm_and_si128(v, mask_0f);
        __m128i hi = _mm_and_si128(_mm_srli_epi16(v, 4), mask_0f);
        __m128i m_lo = _mm_shuffle_epi8(lut_lo, lo);
        __m128i m_hi = _mm_shuffle_epi8(lut_hi, hi);
        __m128i hit  = _mm_and_si128(m_lo, m_hi);
        if (!_mm_testz_si128(hit, hit)) return true;
    }
    return has_degenerate_base_scalar(data + i, n - i);
}

__attribute__((target("avx2")))
static bool has_degenerate_base_avx2(const char* data,
                                     std::size_t n) noexcept {
    // Broadcast 16-byte LUT into both 128-bit lanes (vpshufb is per-lane).
    const __m256i lut_lo = _mm256_broadcastsi128_si256(
        _mm_load_si128(reinterpret_cast<const __m128i*>(kLutLoBitset)));
    const __m256i lut_hi = _mm256_broadcastsi128_si256(
        _mm_load_si128(reinterpret_cast<const __m128i*>(kLutHiRowmask)));
    const __m256i mask_0f = _mm256_set1_epi8(0x0F);
    std::size_t i = 0;
    for (; i + 32 <= n; i += 32) {
        __m256i v  = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(data + i));
        __m256i lo = _mm256_and_si256(v, mask_0f);
        __m256i hi = _mm256_and_si256(_mm256_srli_epi16(v, 4), mask_0f);
        __m256i m_lo = _mm256_shuffle_epi8(lut_lo, lo);
        __m256i m_hi = _mm256_shuffle_epi8(lut_hi, hi);
        __m256i hit  = _mm256_and_si256(m_lo, m_hi);
        if (!_mm256_testz_si256(hit, hit)) return true;
    }
    return has_degenerate_base_scalar(data + i, n - i);
}

__attribute__((target("avx512bw,avx512f")))
static bool has_degenerate_base_avx512bw(const char* data,
                                         std::size_t n) noexcept {
    alignas(64) static const std::int8_t lo_bcast[64] = {
        kLutLoBitset[0],  kLutLoBitset[1],  kLutLoBitset[2],  kLutLoBitset[3],
        kLutLoBitset[4],  kLutLoBitset[5],  kLutLoBitset[6],  kLutLoBitset[7],
        kLutLoBitset[8],  kLutLoBitset[9],  kLutLoBitset[10], kLutLoBitset[11],
        kLutLoBitset[12], kLutLoBitset[13], kLutLoBitset[14], kLutLoBitset[15],
        kLutLoBitset[0],  kLutLoBitset[1],  kLutLoBitset[2],  kLutLoBitset[3],
        kLutLoBitset[4],  kLutLoBitset[5],  kLutLoBitset[6],  kLutLoBitset[7],
        kLutLoBitset[8],  kLutLoBitset[9],  kLutLoBitset[10], kLutLoBitset[11],
        kLutLoBitset[12], kLutLoBitset[13], kLutLoBitset[14], kLutLoBitset[15],
        kLutLoBitset[0],  kLutLoBitset[1],  kLutLoBitset[2],  kLutLoBitset[3],
        kLutLoBitset[4],  kLutLoBitset[5],  kLutLoBitset[6],  kLutLoBitset[7],
        kLutLoBitset[8],  kLutLoBitset[9],  kLutLoBitset[10], kLutLoBitset[11],
        kLutLoBitset[12], kLutLoBitset[13], kLutLoBitset[14], kLutLoBitset[15],
        kLutLoBitset[0],  kLutLoBitset[1],  kLutLoBitset[2],  kLutLoBitset[3],
        kLutLoBitset[4],  kLutLoBitset[5],  kLutLoBitset[6],  kLutLoBitset[7],
        kLutLoBitset[8],  kLutLoBitset[9],  kLutLoBitset[10], kLutLoBitset[11],
        kLutLoBitset[12], kLutLoBitset[13], kLutLoBitset[14], kLutLoBitset[15],
    };
    alignas(64) static const std::int8_t hi_bcast[64] = {
        kLutHiRowmask[0],  kLutHiRowmask[1],  kLutHiRowmask[2],  kLutHiRowmask[3],
        kLutHiRowmask[4],  kLutHiRowmask[5],  kLutHiRowmask[6],  kLutHiRowmask[7],
        kLutHiRowmask[8],  kLutHiRowmask[9],  kLutHiRowmask[10], kLutHiRowmask[11],
        kLutHiRowmask[12], kLutHiRowmask[13], kLutHiRowmask[14], kLutHiRowmask[15],
        kLutHiRowmask[0],  kLutHiRowmask[1],  kLutHiRowmask[2],  kLutHiRowmask[3],
        kLutHiRowmask[4],  kLutHiRowmask[5],  kLutHiRowmask[6],  kLutHiRowmask[7],
        kLutHiRowmask[8],  kLutHiRowmask[9],  kLutHiRowmask[10], kLutHiRowmask[11],
        kLutHiRowmask[12], kLutHiRowmask[13], kLutHiRowmask[14], kLutHiRowmask[15],
        kLutHiRowmask[0],  kLutHiRowmask[1],  kLutHiRowmask[2],  kLutHiRowmask[3],
        kLutHiRowmask[4],  kLutHiRowmask[5],  kLutHiRowmask[6],  kLutHiRowmask[7],
        kLutHiRowmask[8],  kLutHiRowmask[9],  kLutHiRowmask[10], kLutHiRowmask[11],
        kLutHiRowmask[12], kLutHiRowmask[13], kLutHiRowmask[14], kLutHiRowmask[15],
        kLutHiRowmask[0],  kLutHiRowmask[1],  kLutHiRowmask[2],  kLutHiRowmask[3],
        kLutHiRowmask[4],  kLutHiRowmask[5],  kLutHiRowmask[6],  kLutHiRowmask[7],
        kLutHiRowmask[8],  kLutHiRowmask[9],  kLutHiRowmask[10], kLutHiRowmask[11],
        kLutHiRowmask[12], kLutHiRowmask[13], kLutHiRowmask[14], kLutHiRowmask[15],
    };
    const __m512i lut_lo = _mm512_load_si512(lo_bcast);
    const __m512i lut_hi = _mm512_load_si512(hi_bcast);
    const __m512i mask_0f = _mm512_set1_epi8(0x0F);
    std::size_t i = 0;
    for (; i + 64 <= n; i += 64) {
        __m512i v  = _mm512_loadu_si512(
            reinterpret_cast<const void*>(data + i));
        __m512i lo = _mm512_and_si512(v, mask_0f);
        __m512i hi = _mm512_and_si512(_mm512_srli_epi16(v, 4), mask_0f);
        __m512i m_lo = _mm512_shuffle_epi8(lut_lo, lo);
        __m512i m_hi = _mm512_shuffle_epi8(lut_hi, hi);
        __m512i hit  = _mm512_and_si512(m_lo, m_hi);
        __mmask64 m  = _mm512_test_epi8_mask(hit, hit);
        if (m) return true;
    }
    return has_degenerate_base_scalar(data + i, n - i);
}

__attribute__((target("avx512vbmi2,avx512vbmi,avx512bw,avx512f")))
static bool has_degenerate_base_avx512vbmi2(const char* data,
                                            std::size_t n) noexcept {
    // VBMI2 introduces vpcompressb / vpexpandb / vpshrdv. None of those
    // simplifies this kernel further; body matches BW. Separate symbol
    // for tier-effect benchmarking.
    return has_degenerate_base_avx512bw(data, n);
}

#endif // IKAFSSN_DS_X86

// aarch64 kernel (NEON only; SVE/SVE2 fall back to NEON via the dispatcher).
#if IKAFSSN_DS_ARM

static bool has_degenerate_base_neon(const char* data,
                                     std::size_t n) noexcept {
    const uint8x16_t lut_lo = vld1q_u8(
        reinterpret_cast<const std::uint8_t*>(kLutLoBitset));
    const uint8x16_t lut_hi = vld1q_u8(
        reinterpret_cast<const std::uint8_t*>(kLutHiRowmask));
    const uint8x16_t mask_0f = vdupq_n_u8(0x0F);
    std::size_t i = 0;
    for (; i + 16 <= n; i += 16) {
        uint8x16_t v  = vld1q_u8(
            reinterpret_cast<const std::uint8_t*>(data + i));
        uint8x16_t lo = vandq_u8(v, mask_0f);
        uint8x16_t hi = vandq_u8(vshrq_n_u8(v, 4), mask_0f);
        uint8x16_t m_lo = vqtbl1q_u8(lut_lo, lo);
        uint8x16_t m_hi = vqtbl1q_u8(lut_hi, hi);
        uint8x16_t hit  = vandq_u8(m_lo, m_hi);
        if (vmaxvq_u8(hit) != 0) return true;
    }
    return has_degenerate_base_scalar(data + i, n - i);
}

#endif // IKAFSSN_DS_ARM

// Public dispatcher.
bool has_degenerate_base(const char* data, std::size_t n) noexcept {
    if (n == 0) return false;
    SimdCap cap = current_simd_cap();

#if IKAFSSN_DS_X86
    if (cap >= SimdCap::AVX512VBMI2 && n >= kSimdMinBytes_AVX512VBMI2) {
        return has_degenerate_base_avx512vbmi2(data, n);
    }
    if (cap >= SimdCap::AVX512BW && n >= kSimdMinBytes_AVX512BW) {
        return has_degenerate_base_avx512bw(data, n);
    }
    if (cap >= SimdCap::AVX2 && n >= kSimdMinBytes_AVX2) {
        return has_degenerate_base_avx2(data, n);
    }
    if (cap >= SimdCap::SSE42 && n >= kSimdMinBytes_SSE42) {
        return has_degenerate_base_sse42(data, n);
    }
#elif IKAFSSN_DS_ARM
    if (cap >= SimdCap::NEON && n >= kSimdMinBytes_NEON) {
        return has_degenerate_base_neon(data, n);
    }
#else
    (void)cap;
#endif
    return has_degenerate_base_scalar(data, n);
}

} // namespace ikafssn
