#include "core/ncbi2na_unpack.hpp"
#include "util/simd_dispatch.hpp"

#include <cstddef>
#include <cstdint>
#include <cstring>

#if defined(__x86_64__) || defined(__i386__)
  #include <immintrin.h>
  #define IKAFSSN_NU_X86 1
#endif

#if defined(__aarch64__)
  #include <arm_neon.h>
  #define IKAFSSN_NU_ARM 1
  #if __has_include(<arm_sve.h>)
    #include <arm_sve.h>
    #define IKAFSSN_NU_HAVE_SVE_HEADER 1
  #endif
#endif

namespace ikafssn {

// ===========================================================================
// Scalar reference (always compiled): single source of truth and tail handler.
// ===========================================================================
static inline void unpack_scalar(const char*    packed,
                                 std::uint32_t  start,
                                 std::uint32_t  count,
                                 std::uint8_t*  out) noexcept {
    for (std::uint32_t i = 0; i < count; ++i) {
        std::uint32_t pos = start + i;
        std::uint8_t  b   = static_cast<std::uint8_t>(packed[pos >> 2]);
        out[i] = static_cast<std::uint8_t>(
            (b >> (6 - 2 * (pos & 3))) & 0x03);
    }
}

// ===========================================================================
// All x86 / aarch64 SIMD kernels share the same algorithm:
//   1. Load N packed bytes (each packing 4 bases).
//   2. Broadcast each byte to 4 output positions
//      (vpshufb / vpermb / vqtbl1q_u8 / svtbl).
//   3. Apply 3 right shifts (>>6, >>4, >>2) to the broadcast vector.
//      Within each 4-byte group, position j needs the value with that byte
//      right-shifted by (6 - 2*j) before masking with 0x03.
//   4. Pick the correct shift result per-lane via per-position 0x03 masks
//      (mask0 selects pos 0 lanes, mask1 -> pos 1, etc.).
//   5. OR all four contributions and store.
//
// All kernels also assume start_base_pos % 4 == 0 (head-misaligned region is
// sent to the scalar fallback by the dispatcher).
// ===========================================================================

#if IKAFSSN_NU_X86

__attribute__((target("sse4.2,ssse3")))
static void unpack_sse42(const char*    packed,
                         std::uint32_t  start,
                         std::uint32_t  count,
                         std::uint8_t*  out) noexcept {
    const __m128i shuf  = _mm_setr_epi8(0,0,0,0, 1,1,1,1, 2,2,2,2, 3,3,3,3);
    const __m128i mask0 = _mm_set1_epi32(0x00000003);
    const __m128i mask1 = _mm_set1_epi32(0x00000300);
    const __m128i mask2 = _mm_set1_epi32(0x00030000);
    const __m128i mask3 = _mm_set1_epi32(0x03000000);

    std::uint32_t in_byte = start >> 2;
    std::uint32_t i = 0;
    for (; i + 16 <= count; i += 16, in_byte += 4) {
        std::uint32_t in32;
        std::memcpy(&in32, packed + in_byte, 4);
        __m128i in128 = _mm_cvtsi32_si128(static_cast<int>(in32));
        __m128i bcast = _mm_shuffle_epi8(in128, shuf);
        __m128i v6 = _mm_srli_epi16(bcast, 6);
        __m128i v4 = _mm_srli_epi16(bcast, 4);
        __m128i v2 = _mm_srli_epi16(bcast, 2);
        __m128i out128 = _mm_or_si128(
            _mm_or_si128(_mm_and_si128(v6, mask0),
                         _mm_and_si128(v4, mask1)),
            _mm_or_si128(_mm_and_si128(v2, mask2),
                         _mm_and_si128(bcast, mask3)));
        _mm_storeu_si128(reinterpret_cast<__m128i*>(out + i), out128);
    }
    if (i < count) unpack_scalar(packed, start + i, count - i, out + i);
}

__attribute__((target("avx2")))
static void unpack_avx2(const char*    packed,
                        std::uint32_t  start,
                        std::uint32_t  count,
                        std::uint8_t*  out) noexcept {
    // pshufb is per-128-bit-lane; both halves use the same broadcast pattern,
    // and we feed each half with a different 4-byte slice via setr_m128i.
    const __m256i shuf = _mm256_setr_epi8(
        0,0,0,0, 1,1,1,1, 2,2,2,2, 3,3,3,3,
        0,0,0,0, 1,1,1,1, 2,2,2,2, 3,3,3,3);
    const __m256i mask0 = _mm256_set1_epi32(0x00000003);
    const __m256i mask1 = _mm256_set1_epi32(0x00000300);
    const __m256i mask2 = _mm256_set1_epi32(0x00030000);
    const __m256i mask3 = _mm256_set1_epi32(0x03000000);

    std::uint32_t in_byte = start >> 2;
    std::uint32_t i = 0;
    for (; i + 32 <= count; i += 32, in_byte += 8) {
        std::uint64_t in64;
        std::memcpy(&in64, packed + in_byte, 8);
        __m128i lo = _mm_cvtsi32_si128(static_cast<int>(in64 & 0xFFFFFFFFu));
        __m128i hi = _mm_cvtsi32_si128(static_cast<int>(in64 >> 32));
        __m256i in256 = _mm256_setr_m128i(lo, hi);
        __m256i bcast = _mm256_shuffle_epi8(in256, shuf);
        __m256i v6 = _mm256_srli_epi16(bcast, 6);
        __m256i v4 = _mm256_srli_epi16(bcast, 4);
        __m256i v2 = _mm256_srli_epi16(bcast, 2);
        __m256i out256 = _mm256_or_si256(
            _mm256_or_si256(_mm256_and_si256(v6, mask0),
                            _mm256_and_si256(v4, mask1)),
            _mm256_or_si256(_mm256_and_si256(v2, mask2),
                            _mm256_and_si256(bcast, mask3)));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(out + i), out256);
    }
    if (i < count) unpack_scalar(packed, start + i, count - i, out + i);
}

// Common 64-byte broadcast index/shuffle pattern: out[j] = in[j/4].
alignas(64) static constexpr std::uint8_t kBroadcastIdx[64] = {
     0, 0, 0, 0,  1, 1, 1, 1,  2, 2, 2, 2,  3, 3, 3, 3,
     4, 4, 4, 4,  5, 5, 5, 5,  6, 6, 6, 6,  7, 7, 7, 7,
     8, 8, 8, 8,  9, 9, 9, 9, 10,10,10,10, 11,11,11,11,
    12,12,12,12, 13,13,13,13, 14,14,14,14, 15,15,15,15
};

__attribute__((target("avx512bw,avx512f")))
static void unpack_avx512bw(const char*    packed,
                            std::uint32_t  start,
                            std::uint32_t  count,
                            std::uint8_t*  out) noexcept {
    // 16 packed bytes -> broadcast (each byte 4x) -> shift+mask -> 64 bases.
    // vpshufb is per-128-bit-lane: with vbroadcast_i32x4 placing the same
    // 16 bytes into all four 128-bit lanes, each lane's shuf indices select
    // a different 4-byte window.
    const __m512i shuf  = _mm512_load_si512(kBroadcastIdx);
    const __m512i mask0 = _mm512_set1_epi32(0x00000003);
    const __m512i mask1 = _mm512_set1_epi32(0x00000300);
    const __m512i mask2 = _mm512_set1_epi32(0x00030000);
    const __m512i mask3 = _mm512_set1_epi32(0x03000000);

    std::uint32_t in_byte = start >> 2;
    std::uint32_t i = 0;
    for (; i + 64 <= count; i += 64, in_byte += 16) {
        __m128i in128 = _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(packed + in_byte));
        __m512i in512 = _mm512_broadcast_i32x4(in128);
        __m512i bcast = _mm512_shuffle_epi8(in512, shuf);
        __m512i v6 = _mm512_srli_epi16(bcast, 6);
        __m512i v4 = _mm512_srli_epi16(bcast, 4);
        __m512i v2 = _mm512_srli_epi16(bcast, 2);
        __m512i out512 = _mm512_or_si512(
            _mm512_or_si512(_mm512_and_si512(v6, mask0),
                            _mm512_and_si512(v4, mask1)),
            _mm512_or_si512(_mm512_and_si512(v2, mask2),
                            _mm512_and_si512(bcast, mask3)));
        _mm512_storeu_si512(reinterpret_cast<void*>(out + i), out512);
    }
    if (i < count) unpack_scalar(packed, start + i, count - i, out + i);
}

__attribute__((target("avx512vbmi2,avx512vbmi,avx512bw,avx512f")))
static void unpack_avx512vbmi2(const char*    packed,
                               std::uint32_t  start,
                               std::uint32_t  count,
                               std::uint8_t*  out) noexcept {
    // Uses VBMI's vpermb (lane-crossing byte broadcast) under a wider
    // target attribute. VBMI2 introduces vpcompressb/vpexpandb/vpshrdv but
    // no per-byte variable shift, so the body matches the BW shape — the
    // kernel is kept as a separate symbol so the runtime dispatcher can
    // pick the VBMI2 tier without an explicit VBMI ladder.
    const __m512i idx   = _mm512_load_si512(kBroadcastIdx);
    const __m512i mask0 = _mm512_set1_epi32(0x00000003);
    const __m512i mask1 = _mm512_set1_epi32(0x00000300);
    const __m512i mask2 = _mm512_set1_epi32(0x00030000);
    const __m512i mask3 = _mm512_set1_epi32(0x03000000);

    std::uint32_t in_byte = start >> 2;
    std::uint32_t i = 0;
    for (; i + 64 <= count; i += 64, in_byte += 16) {
        __m128i in128 = _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(packed + in_byte));
        __m512i in512 = _mm512_castsi128_si512(in128);
        __m512i bcast = _mm512_permutexvar_epi8(idx, in512);
        __m512i v6 = _mm512_srli_epi16(bcast, 6);
        __m512i v4 = _mm512_srli_epi16(bcast, 4);
        __m512i v2 = _mm512_srli_epi16(bcast, 2);
        __m512i out512 = _mm512_or_si512(
            _mm512_or_si512(_mm512_and_si512(v6, mask0),
                            _mm512_and_si512(v4, mask1)),
            _mm512_or_si512(_mm512_and_si512(v2, mask2),
                            _mm512_and_si512(bcast, mask3)));
        _mm512_storeu_si512(reinterpret_cast<void*>(out + i), out512);
    }
    if (i < count) unpack_scalar(packed, start + i, count - i, out + i);
}

#endif // IKAFSSN_NU_X86

// ===========================================================================
// aarch64.
// ===========================================================================
#if IKAFSSN_NU_ARM

static void unpack_neon(const char*    packed,
                        std::uint32_t  start,
                        std::uint32_t  count,
                        std::uint8_t*  out) noexcept {
    const uint8x16_t shuf  = {0,0,0,0, 1,1,1,1, 2,2,2,2, 3,3,3,3};
    const uint8x16_t mask0 = {0x03,0,0,0, 0x03,0,0,0, 0x03,0,0,0, 0x03,0,0,0};
    const uint8x16_t mask1 = {0,0x03,0,0, 0,0x03,0,0, 0,0x03,0,0, 0,0x03,0,0};
    const uint8x16_t mask2 = {0,0,0x03,0, 0,0,0x03,0, 0,0,0x03,0, 0,0,0x03,0};
    const uint8x16_t mask3 = {0,0,0,0x03, 0,0,0,0x03, 0,0,0,0x03, 0,0,0,0x03};

    std::uint32_t in_byte = start >> 2;
    std::uint32_t i = 0;
    for (; i + 16 <= count; i += 16, in_byte += 4) {
        std::uint32_t in32;
        std::memcpy(&in32, packed + in_byte, 4);
        // Place 4 bytes into lane 0 of a 16-byte register; rest is don't-care.
        uint8x16_t in_vec = vreinterpretq_u8_u32(
            vsetq_lane_u32(in32, vdupq_n_u32(0), 0));
        uint8x16_t bcast = vqtbl1q_u8(in_vec, shuf);
        // NEON has byte-granular immediate right shift.
        uint8x16_t v6 = vshrq_n_u8(bcast, 6);
        uint8x16_t v4 = vshrq_n_u8(bcast, 4);
        uint8x16_t v2 = vshrq_n_u8(bcast, 2);
        uint8x16_t outv = vorrq_u8(
            vorrq_u8(vandq_u8(v6, mask0), vandq_u8(v4, mask1)),
            vorrq_u8(vandq_u8(v2, mask2), vandq_u8(bcast, mask3)));
        vst1q_u8(out + i, outv);
    }
    if (i < count) unpack_scalar(packed, start + i, count - i, out + i);
}

#if IKAFSSN_NU_HAVE_SVE_HEADER

__attribute__((target("+sve")))
static void unpack_sve(const char*    packed,
                       std::uint32_t  start,
                       std::uint32_t  count,
                       std::uint8_t*  out) noexcept {
    // VL-agnostic: process svcntb() bases per iteration. Each iteration
    // consumes svcntb()/4 input bytes (start_base_pos % 4 == 0 by contract).
    const svbool_t  all   = svptrue_b8();
    const svuint8_t lane  = svindex_u8(0, 1);
    const svuint8_t tbl   = svlsr_n_u8_x(all, lane, 2);          // j/4
    const svuint8_t mod4  = svand_n_u8_x(all, lane, 0x03);       // j%4
    const svuint8_t shift = svsub_u8_x(
        all, svdup_n_u8(6),
        svlsl_n_u8_x(all, mod4, 1));                             // 6-2*(j%4)

    const std::uint64_t VL = svcntb();
    std::uint32_t in_byte = start >> 2;
    std::uint32_t i = 0;
    while (i + VL <= count) {
        // Load only the first VL/4 input bytes; trailing lanes are zeroed and
        // never referenced by tbl since tbl[j] = j/4 < VL/4 for all active j.
        svbool_t pg_in = svwhilelt_b8(static_cast<std::uint64_t>(0), VL / 4);
        svuint8_t in_vec = svld1_u8(
            pg_in, reinterpret_cast<const std::uint8_t*>(packed + in_byte));
        svuint8_t bcast = svtbl_u8(in_vec, tbl);
        svuint8_t shifted = svlsr_u8_x(all, bcast, shift);
        svuint8_t result  = svand_n_u8_x(all, shifted, 0x03);
        svst1_u8(all, out + i, result);
        i       += static_cast<std::uint32_t>(VL);
        in_byte += static_cast<std::uint32_t>(VL / 4);
    }
    if (i < count) unpack_scalar(packed, start + i, count - i, out + i);
}

__attribute__((target("+sve2")))
static void unpack_sve2(const char*    packed,
                        std::uint32_t  start,
                        std::uint32_t  count,
                        std::uint8_t*  out) noexcept {
    // SVE2 does not introduce instructions that simplify this kernel further;
    // body matches SVE for tier-effect measurement.
    const svbool_t  all   = svptrue_b8();
    const svuint8_t lane  = svindex_u8(0, 1);
    const svuint8_t tbl   = svlsr_n_u8_x(all, lane, 2);
    const svuint8_t mod4  = svand_n_u8_x(all, lane, 0x03);
    const svuint8_t shift = svsub_u8_x(
        all, svdup_n_u8(6),
        svlsl_n_u8_x(all, mod4, 1));

    const std::uint64_t VL = svcntb();
    std::uint32_t in_byte = start >> 2;
    std::uint32_t i = 0;
    while (i + VL <= count) {
        svbool_t pg_in = svwhilelt_b8(static_cast<std::uint64_t>(0), VL / 4);
        svuint8_t in_vec = svld1_u8(
            pg_in, reinterpret_cast<const std::uint8_t*>(packed + in_byte));
        svuint8_t bcast = svtbl_u8(in_vec, tbl);
        svuint8_t shifted = svlsr_u8_x(all, bcast, shift);
        svuint8_t result  = svand_n_u8_x(all, shifted, 0x03);
        svst1_u8(all, out + i, result);
        i       += static_cast<std::uint32_t>(VL);
        in_byte += static_cast<std::uint32_t>(VL / 4);
    }
    if (i < count) unpack_scalar(packed, start + i, count - i, out + i);
}

#endif // IKAFSSN_NU_HAVE_SVE_HEADER

#endif // IKAFSSN_NU_ARM

// ===========================================================================
// Public dispatcher.
// ===========================================================================
void unpack_ncbi2na_chunk(const char*    packed,
                          std::uint32_t  start,
                          std::uint32_t  count,
                          std::uint8_t*  out) noexcept {
    if (count == 0) return;

    // Send the head-misaligned region to the scalar fallback so the SIMD
    // kernels can assume start % 4 == 0.
    std::uint32_t head = (4u - (start & 3u)) & 3u;
    if (head > count) head = count;
    if (head > 0) {
        unpack_scalar(packed, start, head, out);
        out   += head;
        start += head;
        count -= head;
        if (count == 0) return;
    }

    SimdCap cap = current_simd_cap();
#if IKAFSSN_NU_X86
    if (cap >= SimdCap::AVX512VBMI2 && count >= 64) {
        unpack_avx512vbmi2(packed, start, count, out);
        return;
    }
    if (cap >= SimdCap::AVX512BW && count >= 64) {
        unpack_avx512bw(packed, start, count, out);
        return;
    }
    if (cap >= SimdCap::AVX2 && count >= 32) {
        unpack_avx2(packed, start, count, out);
        return;
    }
    if (cap >= SimdCap::SSE42 && count >= 16) {
        unpack_sse42(packed, start, count, out);
        return;
    }
#elif IKAFSSN_NU_ARM
  #if IKAFSSN_NU_HAVE_SVE_HEADER
    if (cap >= SimdCap::SVE2 && count >= 16) {
        unpack_sve2(packed, start, count, out);
        return;
    }
    if (cap >= SimdCap::SVE && count >= 16) {
        unpack_sve(packed, start, count, out);
        return;
    }
  #endif
    if (cap >= SimdCap::NEON && count >= 16) {
        unpack_neon(packed, start, count, out);
        return;
    }
#else
    (void)cap;
#endif
    unpack_scalar(packed, start, count, out);
}

} // namespace ikafssn
