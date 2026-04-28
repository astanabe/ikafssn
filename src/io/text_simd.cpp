#include "io/text_simd.hpp"
#include "util/simd_dispatch.hpp"

#include <cstddef>
#include <cstdint>

#if defined(__x86_64__) || defined(__i386__)
  #include <immintrin.h>
  #define IKAFSSN_TS_X86 1
#endif

#if defined(__aarch64__)
  #include <arm_neon.h>
  #define IKAFSSN_TS_ARM 1
  #if __has_include(<arm_sve.h>)
    #include <arm_sve.h>
    #define IKAFSSN_TS_HAVE_SVE_HEADER 1
  #endif
#endif

namespace ikafssn {

// ===========================================================================
// Scalar reference (always compiled): handles tail and all fall-back paths.
// ===========================================================================
static inline void toupper_scalar(std::uint8_t* p, std::size_t n) noexcept {
    for (std::size_t i = 0; i < n; ++i) {
        std::uint8_t c = p[i];
        if (c >= 'a' && c <= 'z') p[i] = static_cast<std::uint8_t>(c - 0x20);
    }
}

// ===========================================================================
// x86_64 kernels.
// ===========================================================================
#if IKAFSSN_TS_X86

// SSE4.2 — 16-byte chunks. SSE only has signed compares; use the bias
// (-0x80) trick: subtracting 0x80 maps {'a'..'z'} = {0x61..0x7A} into
// {0xE1..0xFA} as signed bytes ({-31..-6}), so cmpgt against (signed) 0x60-0x80
// = -32 isolates the lowercase range; subtract again against (signed) 0x7A-0x80
// = -6 to bound from above.
//
// Equivalent (and clearer) approach used here: form mask = ((c >= 'a') &
// (c <= 'z')) by combining two greater-than comparisons; then sub 0x20 where
// mask is set.
__attribute__((target("sse4.2")))
static void toupper_sse42(std::uint8_t* p, std::size_t n) noexcept {
    const __m128i va_minus1 = _mm_set1_epi8(static_cast<char>('a' - 1));
    const __m128i vz_plus1  = _mm_set1_epi8(static_cast<char>('z' + 1));
    const __m128i v20       = _mm_set1_epi8(static_cast<char>(0x20));
    const __m128i bias      = _mm_set1_epi8(static_cast<char>(0x80));
    std::size_t i = 0;
    for (; i + 16 <= n; i += 16) {
        __m128i v = _mm_loadu_si128(reinterpret_cast<const __m128i*>(p + i));
        // unsigned compares via signed bias: cmpgt(v - 0x80, 'a'-1 - 0x80)
        __m128i vb = _mm_sub_epi8(v, bias);
        __m128i ge_a = _mm_cmpgt_epi8(vb, _mm_sub_epi8(va_minus1, bias));
        __m128i lt_z1 = _mm_cmpgt_epi8(_mm_sub_epi8(vz_plus1, bias), vb);
        __m128i mask = _mm_and_si128(ge_a, lt_z1);
        __m128i sub  = _mm_and_si128(mask, v20);
        _mm_storeu_si128(reinterpret_cast<__m128i*>(p + i),
                         _mm_sub_epi8(v, sub));
    }
    if (i < n) toupper_scalar(p + i, n - i);
}

// AVX2 — 32-byte chunks. Same shape as SSE4.2 but on __m256i.
__attribute__((target("avx2")))
static void toupper_avx2(std::uint8_t* p, std::size_t n) noexcept {
    const __m256i va_minus1 = _mm256_set1_epi8(static_cast<char>('a' - 1));
    const __m256i vz_plus1  = _mm256_set1_epi8(static_cast<char>('z' + 1));
    const __m256i v20       = _mm256_set1_epi8(static_cast<char>(0x20));
    const __m256i bias      = _mm256_set1_epi8(static_cast<char>(0x80));
    std::size_t i = 0;
    for (; i + 32 <= n; i += 32) {
        __m256i v = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(p + i));
        __m256i vb    = _mm256_sub_epi8(v, bias);
        __m256i ge_a  = _mm256_cmpgt_epi8(vb, _mm256_sub_epi8(va_minus1, bias));
        __m256i lt_z1 = _mm256_cmpgt_epi8(_mm256_sub_epi8(vz_plus1, bias), vb);
        __m256i mask  = _mm256_and_si256(ge_a, lt_z1);
        __m256i sub   = _mm256_and_si256(mask, v20);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(p + i),
                            _mm256_sub_epi8(v, sub));
    }
    if (i < n) toupper_scalar(p + i, n - i);
}

// AVX-512 BW — 64-byte chunks via __mmask64 + masked subtract.
__attribute__((target("avx512bw,avx512f")))
static void toupper_avx512bw(std::uint8_t* p, std::size_t n) noexcept {
    const __m512i va  = _mm512_set1_epi8(static_cast<char>('a'));
    const __m512i vz  = _mm512_set1_epi8(static_cast<char>('z'));
    const __m512i v20 = _mm512_set1_epi8(static_cast<char>(0x20));
    std::size_t i = 0;
    for (; i + 64 <= n; i += 64) {
        __m512i v = _mm512_loadu_si512(reinterpret_cast<const void*>(p + i));
        __mmask64 m_ge = _mm512_cmpge_epu8_mask(v, va);
        __mmask64 m_le = _mm512_cmple_epu8_mask(v, vz);
        __mmask64 m    = m_ge & m_le;
        __m512i out = _mm512_mask_sub_epi8(v, m, v, v20);
        _mm512_storeu_si512(reinterpret_cast<void*>(p + i), out);
    }
    if (i < n) toupper_scalar(p + i, n - i);
}

// AVX-512 VBMI2 — same body, separate symbol for tier-effect measurement.
__attribute__((target("avx512vbmi2,avx512vbmi,avx512bw,avx512f")))
static void toupper_avx512vbmi2(std::uint8_t* p, std::size_t n) noexcept {
    const __m512i va  = _mm512_set1_epi8(static_cast<char>('a'));
    const __m512i vz  = _mm512_set1_epi8(static_cast<char>('z'));
    const __m512i v20 = _mm512_set1_epi8(static_cast<char>(0x20));
    std::size_t i = 0;
    for (; i + 64 <= n; i += 64) {
        __m512i v = _mm512_loadu_si512(reinterpret_cast<const void*>(p + i));
        __mmask64 m_ge = _mm512_cmpge_epu8_mask(v, va);
        __mmask64 m_le = _mm512_cmple_epu8_mask(v, vz);
        __mmask64 m    = m_ge & m_le;
        __m512i out = _mm512_mask_sub_epi8(v, m, v, v20);
        _mm512_storeu_si512(reinterpret_cast<void*>(p + i), out);
    }
    if (i < n) toupper_scalar(p + i, n - i);
}

#endif // IKAFSSN_TS_X86

// ===========================================================================
// aarch64 kernels.
// ===========================================================================
#if IKAFSSN_TS_ARM

// NEON — 16-byte chunks.
static void toupper_neon(std::uint8_t* p, std::size_t n) noexcept {
    const uint8x16_t va  = vdupq_n_u8('a');
    const uint8x16_t vz  = vdupq_n_u8('z');
    const uint8x16_t v20 = vdupq_n_u8(0x20);
    std::size_t i = 0;
    for (; i + 16 <= n; i += 16) {
        uint8x16_t v   = vld1q_u8(p + i);
        uint8x16_t ge  = vcgeq_u8(v, va);
        uint8x16_t le  = vcleq_u8(v, vz);
        uint8x16_t m   = vandq_u8(ge, le);
        uint8x16_t sub = vandq_u8(m, v20);
        vst1q_u8(p + i, vsubq_u8(v, sub));
    }
    if (i < n) toupper_scalar(p + i, n - i);
}

#if IKAFSSN_TS_HAVE_SVE_HEADER

__attribute__((target("sve")))
static void toupper_sve(std::uint8_t* p, std::size_t n) noexcept {
    std::size_t i = 0;
    while (i < n) {
        svbool_t pg = svwhilelt_b8(static_cast<std::uint64_t>(i),
                                   static_cast<std::uint64_t>(n));
        svuint8_t v = svld1_u8(pg, p + i);
        svbool_t ge = svcmpge_n_u8(pg, v, 'a');
        svbool_t le = svcmple_n_u8(pg, v, 'z');
        svbool_t m  = svand_b_z(pg, ge, le);
        svuint8_t out = svsub_n_u8_m(m, v, 0x20);
        svst1_u8(pg, p + i, out);
        i += svcntb();
    }
}

__attribute__((target("sve2,sve")))
static void toupper_sve2(std::uint8_t* p, std::size_t n) noexcept {
    std::size_t i = 0;
    while (i < n) {
        svbool_t pg = svwhilelt_b8(static_cast<std::uint64_t>(i),
                                   static_cast<std::uint64_t>(n));
        svuint8_t v = svld1_u8(pg, p + i);
        svbool_t ge = svcmpge_n_u8(pg, v, 'a');
        svbool_t le = svcmple_n_u8(pg, v, 'z');
        svbool_t m  = svand_b_z(pg, ge, le);
        svuint8_t out = svsub_n_u8_m(m, v, 0x20);
        svst1_u8(pg, p + i, out);
        i += svcntb();
    }
}

#endif // IKAFSSN_TS_HAVE_SVE_HEADER

#endif // IKAFSSN_TS_ARM

// ===========================================================================
// Public dispatcher.
// ===========================================================================
void toupper_inplace_ascii(std::uint8_t* p, std::size_t n) noexcept {
    if (n == 0) return;
    SimdCap cap = current_simd_cap();

#if IKAFSSN_TS_X86
    if (cap >= SimdCap::AVX512VBMI2 && n >= kSimdMinBytes_AVX512VBMI2) {
        toupper_avx512vbmi2(p, n);
        return;
    }
    if (cap >= SimdCap::AVX512BW && n >= kSimdMinBytes_AVX512BW) {
        toupper_avx512bw(p, n);
        return;
    }
    if (cap >= SimdCap::AVX2 && n >= kSimdMinBytes_AVX2) {
        toupper_avx2(p, n);
        return;
    }
    if (cap >= SimdCap::SSE42 && n >= kSimdMinBytes_SSE42) {
        toupper_sse42(p, n);
        return;
    }
#elif IKAFSSN_TS_ARM
  #if IKAFSSN_TS_HAVE_SVE_HEADER
    if (cap >= SimdCap::SVE2 && n >= kSimdMinBytes_SVE2) {
        toupper_sve2(p, n);
        return;
    }
    if (cap >= SimdCap::SVE && n >= kSimdMinBytes_SVE) {
        toupper_sve(p, n);
        return;
    }
  #endif
    if (cap >= SimdCap::NEON && n >= kSimdMinBytes_NEON) {
        toupper_neon(p, n);
        return;
    }
#else
    (void)cap;
#endif
    toupper_scalar(p, n);
}

} // namespace ikafssn
