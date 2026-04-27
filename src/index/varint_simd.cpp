#include "index/varint_simd.hpp"
#include "util/simd_dispatch.hpp"
#include "core/varint.hpp"

#include <cstddef>
#include <cstdint>
#include <cstring>

#if defined(__x86_64__) || defined(__i386__)
  #include <immintrin.h>
  #define IKAFSSN_VS_X86 1
#endif

#if defined(__aarch64__)
  #include <arm_neon.h>
  #define IKAFSSN_VS_ARM 1
  #if __has_include(<arm_sve.h>)
    #include <arm_sve.h>
    #define IKAFSSN_VS_HAVE_SVE_HEADER 1
  #endif
#endif

namespace ikafssn {

// ===========================================================================
// Scalar reference (always compiled). Decodes up to max_count varints from
// [in, in_end). Stops early when running out of input mid-varint. Used as the
// golden path for bit-exact tier comparison and for sub-block tail handling
// inside the SIMD kernels.
// ===========================================================================
static int varint_decode_batch_scalar(const std::uint8_t* in,
                                      const std::uint8_t* in_end,
                                      std::uint32_t* out,
                                      std::size_t* consumed,
                                      int max_count) noexcept {
    const std::uint8_t* p = in;
    int i = 0;
    for (; i < max_count; ++i) {
        if (p >= in_end) break;
        std::uint32_t v = 0;
        unsigned shift = 0;
        std::uint8_t b;
        const std::uint8_t* save = p;
        bool ok = false;
        do {
            if (p >= in_end) {
                // mid-varint truncation: roll back partial consumption
                p = save;
                ok = false;
                break;
            }
            b = *p++;
            v |= static_cast<std::uint32_t>(b & 0x7F) << shift;
            shift += 7;
            if (!(b & 0x80)) { ok = true; break; }
        } while (true);
        if (!ok) break;
        out[i] = v;
    }
    *consumed = static_cast<std::size_t>(p - in);
    return i;
}

// ===========================================================================
// Block-decode helper: given a prefetched continuation-bit mask covering up to
// 64 bytes of input, decode up to max_count varints whose total length fits in
// `block_len` bytes. Returns number of varints decoded; advances *consumed.
//
// Used by every SIMD tier as the kernel of the per-block decode. Inlined into
// each per-target-attribute wrapper so the compiler can keep `cmask` in a
// register and use BMI/CTZ instructions enabled by that target attribute.
// ===========================================================================
static inline int decode_block_with_cmask(const std::uint8_t* in,
                                          std::uint64_t cmask,
                                          int block_len,
                                          std::uint32_t* out,
                                          int max_count,
                                          int* consumed_in_block) noexcept {
    // cmask bit i = 1 iff input byte i has the continuation bit set (bit 7).
    // A varint ends at the first byte position `end_pos` where bit end_pos of
    // cmask is 0. We restrict our scan to the block range [0, block_len).
    int produced = 0;
    int p = 0;
    const std::uint64_t valid_mask = (block_len >= 64)
        ? ~std::uint64_t{0}
        : ((std::uint64_t{1} << block_len) - 1);
    cmask &= valid_mask;
    while (produced < max_count && p < block_len) {
        // Stop bytes are positions [p, block_len) where cmask is 0.
        std::uint64_t shifted_valid = valid_mask & (~std::uint64_t{0} << p);
        std::uint64_t stop_mask = (~cmask) & shifted_valid;
        if (stop_mask == 0) break;  // current varint extends past this block
        int end_pos = __builtin_ctzll(stop_mask);
        int len = end_pos - p + 1;
        if (len > 5) {
            // Malformed varint (> 5 bytes for uint32). Stop here; caller will
            // surface the error on the next read attempt.
            break;
        }
        std::uint32_t v = 0;
        unsigned shift = 0;
        for (int j = 0; j < len; ++j) {
            std::uint8_t b = in[p + j];
            v |= static_cast<std::uint32_t>(b & 0x7F) << shift;
            shift += 7;
        }
        out[produced++] = v;
        p += len;
    }
    *consumed_in_block = p;
    return produced;
}

// ===========================================================================
// x86_64 kernels.
// ===========================================================================
#if IKAFSSN_VS_X86

__attribute__((target("sse4.2")))
static int varint_decode_batch_sse42(const std::uint8_t* in,
                                     const std::uint8_t* in_end,
                                     std::uint32_t* out,
                                     std::size_t* consumed,
                                     int max_count) noexcept {
    const std::uint8_t* p = in;
    int produced = 0;
    while (produced < max_count) {
        std::ptrdiff_t avail = in_end - p;
        if (avail < 16) break;
        __m128i v = _mm_loadu_si128(reinterpret_cast<const __m128i*>(p));
        std::uint32_t cmask32 = static_cast<std::uint32_t>(_mm_movemask_epi8(v)) & 0xFFFFu;
        int consumed_in_block = 0;
        int n = decode_block_with_cmask(p, static_cast<std::uint64_t>(cmask32),
                                        16, out + produced,
                                        max_count - produced, &consumed_in_block);
        if (n == 0 && consumed_in_block == 0) break;
        produced += n;
        p        += consumed_in_block;
    }
    // Tail: scalar to drain the remaining < 16 bytes (or the boundary that the
    // SIMD block declined because it would have crossed 16-byte block).
    std::size_t scalar_consumed = 0;
    int n_tail = varint_decode_batch_scalar(p, in_end, out + produced,
                                            &scalar_consumed,
                                            max_count - produced);
    produced += n_tail;
    p        += scalar_consumed;
    *consumed = static_cast<std::size_t>(p - in);
    return produced;
}

__attribute__((target("avx2,bmi2")))
static int varint_decode_batch_avx2(const std::uint8_t* in,
                                    const std::uint8_t* in_end,
                                    std::uint32_t* out,
                                    std::size_t* consumed,
                                    int max_count) noexcept {
    const std::uint8_t* p = in;
    int produced = 0;
    while (produced < max_count) {
        std::ptrdiff_t avail = in_end - p;
        if (avail < 32) break;
        __m256i v = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(p));
        std::uint64_t cmask = static_cast<std::uint64_t>(
            static_cast<std::uint32_t>(_mm256_movemask_epi8(v)));
        int consumed_in_block = 0;
        int n = decode_block_with_cmask(p, cmask, 32, out + produced,
                                        max_count - produced, &consumed_in_block);
        if (n == 0 && consumed_in_block == 0) break;
        produced += n;
        p        += consumed_in_block;
    }
    std::size_t scalar_consumed = 0;
    int n_tail = varint_decode_batch_scalar(p, in_end, out + produced,
                                            &scalar_consumed,
                                            max_count - produced);
    produced += n_tail;
    p        += scalar_consumed;
    *consumed = static_cast<std::size_t>(p - in);
    return produced;
}

__attribute__((target("avx512bw,avx512f,avx2,bmi2")))
static int varint_decode_batch_avx512bw(const std::uint8_t* in,
                                        const std::uint8_t* in_end,
                                        std::uint32_t* out,
                                        std::size_t* consumed,
                                        int max_count) noexcept {
    const std::uint8_t* p = in;
    int produced = 0;
    while (produced < max_count) {
        std::ptrdiff_t avail = in_end - p;
        if (avail < 64) break;
        __m512i v = _mm512_loadu_si512(reinterpret_cast<const void*>(p));
        // _mm512_movepi8_mask returns a __mmask64 where bit i = sign bit of
        // byte i. The varint continuation bit is exactly that.
        std::uint64_t cmask = _mm512_movepi8_mask(v);
        int consumed_in_block = 0;
        int n = decode_block_with_cmask(p, cmask, 64, out + produced,
                                        max_count - produced, &consumed_in_block);
        if (n == 0 && consumed_in_block == 0) break;
        produced += n;
        p        += consumed_in_block;
    }
    std::size_t scalar_consumed = 0;
    int n_tail = varint_decode_batch_scalar(p, in_end, out + produced,
                                            &scalar_consumed,
                                            max_count - produced);
    produced += n_tail;
    p        += scalar_consumed;
    *consumed = static_cast<std::size_t>(p - in);
    return produced;
}

__attribute__((target("avx512vbmi,avx512bw,avx512f,avx2,bmi2")))
static int varint_decode_batch_avx512vbmi(const std::uint8_t* in,
                                          const std::uint8_t* in_end,
                                          std::uint32_t* out,
                                          std::size_t* consumed,
                                          int max_count) noexcept {
    // VBMI does not add a faster route over BW for this kernel (the bottleneck
    // is the scalar reconstruct). Body matches BW exactly, separate symbol so
    // tier-effect benchmarks can isolate it.
    const std::uint8_t* p = in;
    int produced = 0;
    while (produced < max_count) {
        std::ptrdiff_t avail = in_end - p;
        if (avail < 64) break;
        __m512i v = _mm512_loadu_si512(reinterpret_cast<const void*>(p));
        std::uint64_t cmask = _mm512_movepi8_mask(v);
        int consumed_in_block = 0;
        int n = decode_block_with_cmask(p, cmask, 64, out + produced,
                                        max_count - produced, &consumed_in_block);
        if (n == 0 && consumed_in_block == 0) break;
        produced += n;
        p        += consumed_in_block;
    }
    std::size_t scalar_consumed = 0;
    int n_tail = varint_decode_batch_scalar(p, in_end, out + produced,
                                            &scalar_consumed,
                                            max_count - produced);
    produced += n_tail;
    p        += scalar_consumed;
    *consumed = static_cast<std::size_t>(p - in);
    return produced;
}

__attribute__((target("avx512vbmi2,avx512vbmi,avx512bw,avx512f,avx2,bmi2")))
static int varint_decode_batch_avx512vbmi2(const std::uint8_t* in,
                                           const std::uint8_t* in_end,
                                           std::uint32_t* out,
                                           std::size_t* consumed,
                                           int max_count) noexcept {
    // VBMI2 introduces vpcompressb / vpexpandb / vpshrdv. None of those
    // unlocks a strictly faster shape than BW for the present scalar-
    // reconstruct strategy, so the body matches BW. Kept as a separate symbol
    // so tier-effect benchmarks can isolate it.
    const std::uint8_t* p = in;
    int produced = 0;
    while (produced < max_count) {
        std::ptrdiff_t avail = in_end - p;
        if (avail < 64) break;
        __m512i v = _mm512_loadu_si512(reinterpret_cast<const void*>(p));
        std::uint64_t cmask = _mm512_movepi8_mask(v);
        int consumed_in_block = 0;
        int n = decode_block_with_cmask(p, cmask, 64, out + produced,
                                        max_count - produced, &consumed_in_block);
        if (n == 0 && consumed_in_block == 0) break;
        produced += n;
        p        += consumed_in_block;
    }
    std::size_t scalar_consumed = 0;
    int n_tail = varint_decode_batch_scalar(p, in_end, out + produced,
                                            &scalar_consumed,
                                            max_count - produced);
    produced += n_tail;
    p        += scalar_consumed;
    *consumed = static_cast<std::size_t>(p - in);
    return produced;
}

#endif // IKAFSSN_VS_X86

// ===========================================================================
// aarch64 kernels.
// ===========================================================================
#if IKAFSSN_VS_ARM

// NEON does not have a direct movemask. Build a 16-bit cmask by shifting each
// byte's bit 7 into a unique position and horizontally adding via paired
// shift+add.
static inline std::uint16_t neon_cmask16(uint8x16_t v) {
    // bit7 of each byte -> 1 in a 0/1 mask
    uint8x16_t hi = vshrq_n_u8(v, 7);
    // Multiply each byte by 1<<lane to spread into a single 16-bit word.
    static const uint8x16_t kPow2 = {
        1, 2, 4, 8, 16, 32, 64, 128,
        1, 2, 4, 8, 16, 32, 64, 128
    };
    uint8x16_t weighted = vmulq_u8(hi, kPow2);
    // Sum the two 8-byte halves separately to get two bytes, each holding the
    // popcount-bit-pattern of one half.
    std::uint16_t lo_byte = static_cast<std::uint16_t>(vaddv_u8(vget_low_u8(weighted)));
    std::uint16_t hi_byte = static_cast<std::uint16_t>(vaddv_u8(vget_high_u8(weighted)));
    return static_cast<std::uint16_t>(lo_byte | (hi_byte << 8));
}

static int varint_decode_batch_neon(const std::uint8_t* in,
                                    const std::uint8_t* in_end,
                                    std::uint32_t* out,
                                    std::size_t* consumed,
                                    int max_count) noexcept {
    const std::uint8_t* p = in;
    int produced = 0;
    while (produced < max_count) {
        std::ptrdiff_t avail = in_end - p;
        if (avail < 16) break;
        uint8x16_t v = vld1q_u8(p);
        std::uint64_t cmask = neon_cmask16(v);
        int consumed_in_block = 0;
        int n = decode_block_with_cmask(p, cmask, 16, out + produced,
                                        max_count - produced, &consumed_in_block);
        if (n == 0 && consumed_in_block == 0) break;
        produced += n;
        p        += consumed_in_block;
    }
    std::size_t scalar_consumed = 0;
    int n_tail = varint_decode_batch_scalar(p, in_end, out + produced,
                                            &scalar_consumed,
                                            max_count - produced);
    produced += n_tail;
    p        += scalar_consumed;
    *consumed = static_cast<std::size_t>(p - in);
    return produced;
}

#if IKAFSSN_VS_HAVE_SVE_HEADER

__attribute__((target("sve")))
static int varint_decode_batch_sve(const std::uint8_t* in,
                                   const std::uint8_t* in_end,
                                   std::uint32_t* out,
                                   std::size_t* consumed,
                                   int max_count) noexcept {
    // VL-agnostic: process up to min(svcntb(), 64) bytes per block. Cap at 64
    // to keep cmask in a uint64. For VL > 64 we still emit a 64-bit cmask
    // covering the first 64 bytes.
    const std::uint64_t VL = svcntb();
    const std::uint64_t block = (VL < 64) ? VL : 64;
    const svbool_t pg_block = svwhilelt_b8(static_cast<std::uint64_t>(0), block);

    const std::uint8_t* p = in;
    int produced = 0;
    while (produced < max_count) {
        std::ptrdiff_t avail = in_end - p;
        if (avail < static_cast<std::ptrdiff_t>(block)) break;
        svuint8_t v = svld1_u8(pg_block, p);
        // bit 7 of each byte
        svuint8_t hi = svlsr_n_u8_x(pg_block, v, 7);
        svbool_t cont = svcmpne_n_u8(pg_block, hi, 0);
        // svptest doesn't yield a packed 64-bit mask directly; compute byte-by-
        // byte by extracting predicate via svdup + svand.
        // Straightforward path: build cmask scalarly by unpacking.
        std::uint8_t tmp[64];
        svst1_u8(pg_block, tmp, hi);
        std::uint64_t cmask = 0;
        for (std::uint64_t i = 0; i < block; ++i) {
            if (tmp[i]) cmask |= (std::uint64_t{1} << i);
        }
        (void)cont;
        int consumed_in_block = 0;
        int n = decode_block_with_cmask(p, cmask, static_cast<int>(block),
                                        out + produced,
                                        max_count - produced, &consumed_in_block);
        if (n == 0 && consumed_in_block == 0) break;
        produced += n;
        p        += consumed_in_block;
    }
    std::size_t scalar_consumed = 0;
    int n_tail = varint_decode_batch_scalar(p, in_end, out + produced,
                                            &scalar_consumed,
                                            max_count - produced);
    produced += n_tail;
    p        += scalar_consumed;
    *consumed = static_cast<std::size_t>(p - in);
    return produced;
}

__attribute__((target("sve2,sve")))
static int varint_decode_batch_sve2(const std::uint8_t* in,
                                    const std::uint8_t* in_end,
                                    std::uint32_t* out,
                                    std::size_t* consumed,
                                    int max_count) noexcept {
    // SVE2 adds bext but our reconstruct is scalar; body matches SVE.
    const std::uint64_t VL = svcntb();
    const std::uint64_t block = (VL < 64) ? VL : 64;
    const svbool_t pg_block = svwhilelt_b8(static_cast<std::uint64_t>(0), block);

    const std::uint8_t* p = in;
    int produced = 0;
    while (produced < max_count) {
        std::ptrdiff_t avail = in_end - p;
        if (avail < static_cast<std::ptrdiff_t>(block)) break;
        svuint8_t v = svld1_u8(pg_block, p);
        svuint8_t hi = svlsr_n_u8_x(pg_block, v, 7);
        std::uint8_t tmp[64];
        svst1_u8(pg_block, tmp, hi);
        std::uint64_t cmask = 0;
        for (std::uint64_t i = 0; i < block; ++i) {
            if (tmp[i]) cmask |= (std::uint64_t{1} << i);
        }
        int consumed_in_block = 0;
        int n = decode_block_with_cmask(p, cmask, static_cast<int>(block),
                                        out + produced,
                                        max_count - produced, &consumed_in_block);
        if (n == 0 && consumed_in_block == 0) break;
        produced += n;
        p        += consumed_in_block;
    }
    std::size_t scalar_consumed = 0;
    int n_tail = varint_decode_batch_scalar(p, in_end, out + produced,
                                            &scalar_consumed,
                                            max_count - produced);
    produced += n_tail;
    p        += scalar_consumed;
    *consumed = static_cast<std::size_t>(p - in);
    return produced;
}

#endif // IKAFSSN_VS_HAVE_SVE_HEADER

#endif // IKAFSSN_VS_ARM

// ===========================================================================
// Public dispatcher.
// ===========================================================================
int varint_decode_batch(const std::uint8_t* in,
                        const std::uint8_t* in_end,
                        std::uint32_t* out_values,
                        std::size_t* out_consumed,
                        int max_count) noexcept {
    if (in >= in_end || max_count <= 0) {
        *out_consumed = 0;
        return 0;
    }
    if (max_count > kVarintMaxBatch) max_count = kVarintMaxBatch;

    SimdCap cap = current_simd_cap();
#if IKAFSSN_VS_X86
    if (cap >= SimdCap::AVX512VBMI2 && (in_end - in) >= 64) {
        return varint_decode_batch_avx512vbmi2(in, in_end, out_values,
                                               out_consumed, max_count);
    }
    if (cap >= SimdCap::AVX512VBMI && (in_end - in) >= 64) {
        return varint_decode_batch_avx512vbmi(in, in_end, out_values,
                                              out_consumed, max_count);
    }
    if (cap >= SimdCap::AVX512BW && (in_end - in) >= 64) {
        return varint_decode_batch_avx512bw(in, in_end, out_values,
                                            out_consumed, max_count);
    }
    if (cap >= SimdCap::AVX2 && (in_end - in) >= 32 &&
        !current_simd_flags().slow_bmi2) {
        return varint_decode_batch_avx2(in, in_end, out_values,
                                        out_consumed, max_count);
    }
    if (cap >= SimdCap::SSE42 && (in_end - in) >= 16) {
        return varint_decode_batch_sse42(in, in_end, out_values,
                                         out_consumed, max_count);
    }
#elif IKAFSSN_VS_ARM
  #if IKAFSSN_VS_HAVE_SVE_HEADER
    if (cap >= SimdCap::SVE2 && (in_end - in) >= 16) {
        return varint_decode_batch_sve2(in, in_end, out_values,
                                        out_consumed, max_count);
    }
    if (cap >= SimdCap::SVE && (in_end - in) >= 16) {
        return varint_decode_batch_sve(in, in_end, out_values,
                                       out_consumed, max_count);
    }
  #endif
    if (cap >= SimdCap::NEON && (in_end - in) >= 16) {
        return varint_decode_batch_neon(in, in_end, out_values,
                                        out_consumed, max_count);
    }
#else
    (void)cap;
#endif
    return varint_decode_batch_scalar(in, in_end, out_values,
                                      out_consumed, max_count);
}

} // namespace ikafssn
