#include "core/spaced_seed_simd.hpp"
#include "core/spaced_seed.hpp"
#include "util/simd_dispatch.hpp"

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <type_traits>

#if defined(__x86_64__) || defined(__i386__)
  #include <immintrin.h>
  #define IKAFSSN_SS_X86 1
#endif

#if defined(__aarch64__)
  #include <arm_neon.h>
  #define IKAFSSN_SS_ARM 1
#endif

namespace ikafssn {

// ===========================================================================
// Scalar reference: bit-exact source of truth.
// ===========================================================================
template <SpacedSeedExtractor Ext, typename KmerInt>
__attribute__((always_inline)) static inline
void extract_kmer_ct_batch_scalar(const std::uint64_t* a, std::size_t n,
                                  KmerInt* out) noexcept {
    for (std::size_t i = 0; i < n; ++i) {
        out[i] = extract_kmer_ct<Ext, KmerInt>(a[i]);
    }
}

// ===========================================================================
// x86_64 kernels.
//
// All four tiers share the same shape: load a chunk of u64 accumulators,
// apply each contiguous-run mask+shift, accumulate via OR, narrow to KmerInt,
// and store. Run application is unrolled by an `if constexpr`-recursive
// template (one per tier — the recursion plus the per-tier target attribute
// keep every shift count visible to the codegen as a compile-time immediate).
//
// AVX-512BW and AVX-512VBMI2 share the body — the latter is a separate
// symbol for tier-effect benchmarking. VBMI2 introduces vpcompressb /
// vpexpandb / vpshrdv but none of those simplifies this kernel further.
// ===========================================================================
#if IKAFSSN_SS_X86

// SSE4.2 — recursive run application + chunk loop.
template <SpacedSeedExtractor Ext, std::size_t I>
__attribute__((target("sse4.2,ssse3"), always_inline)) static inline
__m128i apply_runs_sse42(__m128i kmer, __m128i accum) noexcept {
    if constexpr (I < Ext.num_runs) {
        constexpr int as = Ext.runs[I].accum_shift;
        constexpr int ks = Ext.runs[I].kmer_shift;
        constexpr std::int64_t m =
            static_cast<std::int64_t>((1ULL << Ext.runs[I].width_bits) - 1);
        __m128i shifted = _mm_srli_epi64(accum, as);
        __m128i masked  = _mm_and_si128(shifted, _mm_set1_epi64x(m));
        __m128i shl     = _mm_slli_epi64(masked, ks);
        kmer = _mm_or_si128(kmer, shl);
        return apply_runs_sse42<Ext, I + 1>(kmer, accum);
    } else {
        return kmer;
    }
}

template <SpacedSeedExtractor Ext, typename KmerInt>
__attribute__((target("sse4.2,ssse3")))
static void extract_kmer_ct_batch_sse42(const std::uint64_t* a, std::size_t n,
                                        KmerInt* out) noexcept {
    std::size_t i = 0;
    for (; i + 2 <= n; i += 2) {
        __m128i v = _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(a + i));
        __m128i kmer = apply_runs_sse42<Ext, 0>(_mm_setzero_si128(), v);
        if constexpr (sizeof(KmerInt) == 4) {
            // [low32(u64_0), high32(u64_0), low32(u64_1), high32(u64_1)]
            // -> shuffle to [low0, low1, low0, low1], store low 64 bits.
            __m128i p = _mm_shuffle_epi32(kmer, _MM_SHUFFLE(2, 0, 2, 0));
            _mm_storel_epi64(reinterpret_cast<__m128i*>(out + i), p);
        } else {
            // u16: scalar narrow (only 2 lanes).
            alignas(16) std::uint64_t buf[2];
            _mm_store_si128(reinterpret_cast<__m128i*>(buf), kmer);
            out[i + 0] = static_cast<std::uint16_t>(buf[0]);
            out[i + 1] = static_cast<std::uint16_t>(buf[1]);
        }
    }
    for (; i < n; ++i) out[i] = extract_kmer_ct<Ext, KmerInt>(a[i]);
}

// AVX2 — recursive run application + chunk loop.
template <SpacedSeedExtractor Ext, std::size_t I>
__attribute__((target("avx2"), always_inline)) static inline
__m256i apply_runs_avx2(__m256i kmer, __m256i accum) noexcept {
    if constexpr (I < Ext.num_runs) {
        constexpr int as = Ext.runs[I].accum_shift;
        constexpr int ks = Ext.runs[I].kmer_shift;
        constexpr std::int64_t m =
            static_cast<std::int64_t>((1ULL << Ext.runs[I].width_bits) - 1);
        __m256i shifted = _mm256_srli_epi64(accum, as);
        __m256i masked  = _mm256_and_si256(shifted, _mm256_set1_epi64x(m));
        __m256i shl     = _mm256_slli_epi64(masked, ks);
        kmer = _mm256_or_si256(kmer, shl);
        return apply_runs_avx2<Ext, I + 1>(kmer, accum);
    } else {
        return kmer;
    }
}

template <SpacedSeedExtractor Ext, typename KmerInt>
__attribute__((target("avx2")))
static void extract_kmer_ct_batch_avx2(const std::uint64_t* a, std::size_t n,
                                       KmerInt* out) noexcept {
    std::size_t i = 0;
    for (; i + 4 <= n; i += 4) {
        __m256i v = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(a + i));
        __m256i kmer = apply_runs_avx2<Ext, 0>(_mm256_setzero_si256(), v);
        if constexpr (sizeof(KmerInt) == 4) {
            // Gather low32 of each u64 (u32 lanes 0,2,4,6) into low 128 bits.
            const __m256i pack_idx =
                _mm256_setr_epi32(0, 2, 4, 6, 0, 0, 0, 0);
            __m256i packed = _mm256_permutevar8x32_epi32(kmer, pack_idx);
            _mm_storeu_si128(reinterpret_cast<__m128i*>(out + i),
                             _mm256_castsi256_si128(packed));
        } else {
            alignas(32) std::uint64_t buf[4];
            _mm256_store_si256(reinterpret_cast<__m256i*>(buf), kmer);
            out[i + 0] = static_cast<std::uint16_t>(buf[0]);
            out[i + 1] = static_cast<std::uint16_t>(buf[1]);
            out[i + 2] = static_cast<std::uint16_t>(buf[2]);
            out[i + 3] = static_cast<std::uint16_t>(buf[3]);
        }
    }
    for (; i < n; ++i) out[i] = extract_kmer_ct<Ext, KmerInt>(a[i]);
}

// AVX-512BW — recursive run application + chunk loop.
template <SpacedSeedExtractor Ext, std::size_t I>
__attribute__((target("avx512bw,avx512f"), always_inline)) static inline
__m512i apply_runs_avx512bw(__m512i kmer, __m512i accum) noexcept {
    if constexpr (I < Ext.num_runs) {
        constexpr int as = Ext.runs[I].accum_shift;
        constexpr int ks = Ext.runs[I].kmer_shift;
        constexpr std::int64_t m =
            static_cast<std::int64_t>((1ULL << Ext.runs[I].width_bits) - 1);
        __m512i shifted = _mm512_srli_epi64(accum, as);
        __m512i masked  = _mm512_and_si512(shifted, _mm512_set1_epi64(m));
        __m512i shl     = _mm512_slli_epi64(masked, ks);
        kmer = _mm512_or_si512(kmer, shl);
        return apply_runs_avx512bw<Ext, I + 1>(kmer, accum);
    } else {
        return kmer;
    }
}

template <SpacedSeedExtractor Ext, typename KmerInt>
__attribute__((target("avx512bw,avx512f")))
static void extract_kmer_ct_batch_avx512bw(const std::uint64_t* a,
                                           std::size_t n,
                                           KmerInt* out) noexcept {
    std::size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        __m512i v = _mm512_loadu_si512(
            reinterpret_cast<const void*>(a + i));
        __m512i kmer = apply_runs_avx512bw<Ext, 0>(_mm512_setzero_si512(), v);
        if constexpr (sizeof(KmerInt) == 4) {
            __m256i out256 = _mm512_cvtepi64_epi32(kmer);
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(out + i), out256);
        } else {
            __m128i out128 = _mm512_cvtepi64_epi16(kmer);
            _mm_storeu_si128(reinterpret_cast<__m128i*>(out + i), out128);
        }
    }
    for (; i < n; ++i) out[i] = extract_kmer_ct<Ext, KmerInt>(a[i]);
}

// AVX-512VBMI2 — separate target/symbol, body matches BW.
template <SpacedSeedExtractor Ext, std::size_t I>
__attribute__((target("avx512vbmi2,avx512vbmi,avx512bw,avx512f"),
               always_inline)) static inline
__m512i apply_runs_avx512vbmi2(__m512i kmer, __m512i accum) noexcept {
    if constexpr (I < Ext.num_runs) {
        constexpr int as = Ext.runs[I].accum_shift;
        constexpr int ks = Ext.runs[I].kmer_shift;
        constexpr std::int64_t m =
            static_cast<std::int64_t>((1ULL << Ext.runs[I].width_bits) - 1);
        __m512i shifted = _mm512_srli_epi64(accum, as);
        __m512i masked  = _mm512_and_si512(shifted, _mm512_set1_epi64(m));
        __m512i shl     = _mm512_slli_epi64(masked, ks);
        kmer = _mm512_or_si512(kmer, shl);
        return apply_runs_avx512vbmi2<Ext, I + 1>(kmer, accum);
    } else {
        return kmer;
    }
}

template <SpacedSeedExtractor Ext, typename KmerInt>
__attribute__((target("avx512vbmi2,avx512vbmi,avx512bw,avx512f")))
static void extract_kmer_ct_batch_avx512vbmi2(const std::uint64_t* a,
                                              std::size_t n,
                                              KmerInt* out) noexcept {
    std::size_t i = 0;
    for (; i + 8 <= n; i += 8) {
        __m512i v = _mm512_loadu_si512(
            reinterpret_cast<const void*>(a + i));
        __m512i kmer = apply_runs_avx512vbmi2<Ext, 0>(
            _mm512_setzero_si512(), v);
        if constexpr (sizeof(KmerInt) == 4) {
            __m256i out256 = _mm512_cvtepi64_epi32(kmer);
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(out + i), out256);
        } else {
            __m128i out128 = _mm512_cvtepi64_epi16(kmer);
            _mm_storeu_si128(reinterpret_cast<__m128i*>(out + i), out128);
        }
    }
    for (; i < n; ++i) out[i] = extract_kmer_ct<Ext, KmerInt>(a[i]);
}

#endif // IKAFSSN_SS_X86

// ===========================================================================
// aarch64 NEON kernel.
//
// Same shape as the x86_64 SSE4.2 path: 2 u64 lanes per iteration in a Q
// register. NEON shift-by-immediate intrinsics (vshrq_n_u64 / vshlq_n_u64)
// require a compile-time constant, which is satisfied because each
// Ext.runs[I].accum_shift / kmer_shift comes from a constexpr NTTP. The
// shift count must additionally be in [1, 63]; an `if constexpr` guards
// the zero-shift case (vshrq_n_u64 rejects 0, and a zero shift is a no-op).
// SVE / SVE2 are not used here — Cortex-A72 has no SVE, and a SVE-native
// extractor would be orthogonal hand-coded work; the NEON tier covers all
// aarch64 CPUs from the ASIMD baseline upward at runtime.
// ===========================================================================
#if IKAFSSN_SS_ARM

template <SpacedSeedExtractor Ext, std::size_t I>
__attribute__((always_inline)) static inline
uint64x2_t apply_runs_neon(uint64x2_t kmer, uint64x2_t accum) noexcept {
    if constexpr (I < Ext.num_runs) {
        constexpr int as = Ext.runs[I].accum_shift;
        constexpr int ks = Ext.runs[I].kmer_shift;
        constexpr std::uint64_t m =
            (1ULL << Ext.runs[I].width_bits) - 1;
        uint64x2_t shifted;
        if constexpr (as == 0) {
            shifted = accum;
        } else {
            shifted = vshrq_n_u64(accum, as);
        }
        uint64x2_t masked = vandq_u64(shifted, vdupq_n_u64(m));
        uint64x2_t shl;
        if constexpr (ks == 0) {
            shl = masked;
        } else {
            shl = vshlq_n_u64(masked, ks);
        }
        kmer = vorrq_u64(kmer, shl);
        return apply_runs_neon<Ext, I + 1>(kmer, accum);
    } else {
        return kmer;
    }
}

template <SpacedSeedExtractor Ext, typename KmerInt>
static void extract_kmer_ct_batch_neon(const std::uint64_t* a, std::size_t n,
                                       KmerInt* out) noexcept {
    std::size_t i = 0;
    for (; i + 2 <= n; i += 2) {
        uint64x2_t v = vld1q_u64(a + i);
        uint64x2_t kmer = apply_runs_neon<Ext, 0>(vdupq_n_u64(0), v);
        if constexpr (sizeof(KmerInt) == 4) {
            // u64x2 -> u32x2 (truncate low 32 bits of each lane).
            uint32x2_t p = vmovn_u64(kmer);
            vst1_u32(reinterpret_cast<std::uint32_t*>(out + i), p);
        } else {
            // u16: scalar narrow (only 2 lanes).
            std::uint64_t buf[2];
            vst1q_u64(buf, kmer);
            out[i + 0] = static_cast<std::uint16_t>(buf[0]);
            out[i + 1] = static_cast<std::uint16_t>(buf[1]);
        }
    }
    for (; i < n; ++i) out[i] = extract_kmer_ct<Ext, KmerInt>(a[i]);
}

#endif // IKAFSSN_SS_ARM

// ===========================================================================
// Per-Ext dispatcher: extract_kmer_ct_batch<Ext, KmerInt>.
//
// 24 explicit specializations are emitted via the X-macro below. Each
// dispatcher selects the highest available tier whose chunk-width fits n,
// then falls back to the always-present scalar path. n thresholds match the
// per-tier chunk widths so the SIMD path is never invoked when its loop body
// would not run a single iteration.
// ===========================================================================

#define IKAFSSN_FOR_EACH_SEED_EXT(M)                                           \
    M(EXT_K8_T13_COD,  std::uint16_t)                                          \
    M(EXT_K8_T13_OPT,  std::uint16_t)                                          \
    M(EXT_K8_T15_COD,  std::uint16_t)                                          \
    M(EXT_K8_T15_OPT,  std::uint16_t)                                          \
    M(EXT_K8_T18_COD,  std::uint16_t)                                          \
    M(EXT_K8_T18_OPT,  std::uint16_t)                                          \
    M(EXT_K9_T13_COD,  std::uint32_t)                                          \
    M(EXT_K9_T13_OPT,  std::uint32_t)                                          \
    M(EXT_K9_T15_COD,  std::uint32_t)                                          \
    M(EXT_K9_T15_OPT,  std::uint32_t)                                          \
    M(EXT_K9_T18_COD,  std::uint32_t)                                          \
    M(EXT_K9_T18_OPT,  std::uint32_t)                                          \
    M(EXT_K11_T16_COD, std::uint32_t)                                          \
    M(EXT_K11_T16_OPT, std::uint32_t)                                          \
    M(EXT_K11_T18_COD, std::uint32_t)                                          \
    M(EXT_K11_T18_OPT, std::uint32_t)                                          \
    M(EXT_K11_T21_COD, std::uint32_t)                                          \
    M(EXT_K11_T21_OPT, std::uint32_t)                                          \
    M(EXT_K12_T16_COD, std::uint32_t)                                          \
    M(EXT_K12_T16_OPT, std::uint32_t)                                          \
    M(EXT_K12_T18_COD, std::uint32_t)                                          \
    M(EXT_K12_T18_OPT, std::uint32_t)                                          \
    M(EXT_K12_T21_COD, std::uint32_t)                                          \
    M(EXT_K12_T21_OPT, std::uint32_t)

#if IKAFSSN_SS_X86
#define IKAFSSN_DEFINE_SEED_DISPATCHER(EXT_VAL, KmerT)                         \
template <>                                                                     \
void extract_kmer_ct_batch<EXT_VAL, KmerT>(const std::uint64_t* a,              \
                                            std::size_t n,                      \
                                            KmerT* out) noexcept {              \
    if (n == 0) return;                                                         \
    SimdCap cap = current_simd_cap();                                           \
    if (cap >= SimdCap::AVX512VBMI2 && n >= 8) {                                \
        extract_kmer_ct_batch_avx512vbmi2<EXT_VAL, KmerT>(a, n, out);           \
        return;                                                                 \
    }                                                                           \
    if (cap >= SimdCap::AVX512BW && n >= 8) {                                   \
        extract_kmer_ct_batch_avx512bw<EXT_VAL, KmerT>(a, n, out);              \
        return;                                                                 \
    }                                                                           \
    if (cap >= SimdCap::AVX2 && n >= 4) {                                       \
        extract_kmer_ct_batch_avx2<EXT_VAL, KmerT>(a, n, out);                  \
        return;                                                                 \
    }                                                                           \
    if (cap >= SimdCap::SSE42 && n >= 2) {                                      \
        extract_kmer_ct_batch_sse42<EXT_VAL, KmerT>(a, n, out);                 \
        return;                                                                 \
    }                                                                           \
    extract_kmer_ct_batch_scalar<EXT_VAL, KmerT>(a, n, out);                    \
}
#elif IKAFSSN_SS_ARM
#define IKAFSSN_DEFINE_SEED_DISPATCHER(EXT_VAL, KmerT)                         \
template <>                                                                     \
void extract_kmer_ct_batch<EXT_VAL, KmerT>(const std::uint64_t* a,              \
                                            std::size_t n,                      \
                                            KmerT* out) noexcept {              \
    if (n == 0) return;                                                         \
    SimdCap cap = current_simd_cap();                                           \
    if (cap >= SimdCap::NEON && n >= 2) {                                       \
        extract_kmer_ct_batch_neon<EXT_VAL, KmerT>(a, n, out);                  \
        return;                                                                 \
    }                                                                           \
    extract_kmer_ct_batch_scalar<EXT_VAL, KmerT>(a, n, out);                    \
}
#else
#define IKAFSSN_DEFINE_SEED_DISPATCHER(EXT_VAL, KmerT)                         \
template <>                                                                     \
void extract_kmer_ct_batch<EXT_VAL, KmerT>(const std::uint64_t* a,              \
                                            std::size_t n,                      \
                                            KmerT* out) noexcept {              \
    if (n == 0) return;                                                         \
    extract_kmer_ct_batch_scalar<EXT_VAL, KmerT>(a, n, out);                    \
}
#endif

IKAFSSN_FOR_EACH_SEED_EXT(IKAFSSN_DEFINE_SEED_DISPATCHER)

#undef IKAFSSN_DEFINE_SEED_DISPATCHER

// ===========================================================================
// Mask-runtime entry point. Switch on the runtime mask to pick the matching
// NTTP specialization. Same shape as extract_for_mask().
// ===========================================================================
template <typename KmerInt>
void extract_for_mask_batch(const std::uint64_t* a, std::size_t n,
                            std::uint32_t mask, KmerInt* out) noexcept {
    if (n == 0) return;

    auto fallback = [&]() {
        // Custom mask not in the canonical 24: fall back to per-element
        // extract_for_mask(), which itself recovers the extractor at runtime.
        for (std::size_t i = 0; i < n; ++i) {
            out[i] = extract_for_mask<KmerInt>(a[i], mask);
        }
    };

    if constexpr (std::is_same_v<KmerInt, std::uint16_t>) {
        switch (mask) {
            case MASK_K8_T13_CODING:
                extract_kmer_ct_batch<EXT_K8_T13_COD, KmerInt>(a, n, out); return;
            case MASK_K8_T13_OPTIMAL:
                extract_kmer_ct_batch<EXT_K8_T13_OPT, KmerInt>(a, n, out); return;
            case MASK_K8_T15_CODING:
                extract_kmer_ct_batch<EXT_K8_T15_COD, KmerInt>(a, n, out); return;
            case MASK_K8_T15_OPTIMAL:
                extract_kmer_ct_batch<EXT_K8_T15_OPT, KmerInt>(a, n, out); return;
            case MASK_K8_T18_CODING:
                extract_kmer_ct_batch<EXT_K8_T18_COD, KmerInt>(a, n, out); return;
            case MASK_K8_T18_OPTIMAL:
                extract_kmer_ct_batch<EXT_K8_T18_OPT, KmerInt>(a, n, out); return;
            default: fallback(); return;
        }
    } else if constexpr (std::is_same_v<KmerInt, std::uint32_t>) {
        switch (mask) {
            case MASK_K9_T13_CODING:
                extract_kmer_ct_batch<EXT_K9_T13_COD, KmerInt>(a, n, out); return;
            case MASK_K9_T13_OPTIMAL:
                extract_kmer_ct_batch<EXT_K9_T13_OPT, KmerInt>(a, n, out); return;
            case MASK_K9_T15_CODING:
                extract_kmer_ct_batch<EXT_K9_T15_COD, KmerInt>(a, n, out); return;
            case MASK_K9_T15_OPTIMAL:
                extract_kmer_ct_batch<EXT_K9_T15_OPT, KmerInt>(a, n, out); return;
            case MASK_K9_T18_CODING:
                extract_kmer_ct_batch<EXT_K9_T18_COD, KmerInt>(a, n, out); return;
            case MASK_K9_T18_OPTIMAL:
                extract_kmer_ct_batch<EXT_K9_T18_OPT, KmerInt>(a, n, out); return;
            case MASK_K11_T16_CODING:
                extract_kmer_ct_batch<EXT_K11_T16_COD, KmerInt>(a, n, out); return;
            case MASK_K11_T16_OPTIMAL:
                extract_kmer_ct_batch<EXT_K11_T16_OPT, KmerInt>(a, n, out); return;
            case MASK_K11_T18_CODING:
                extract_kmer_ct_batch<EXT_K11_T18_COD, KmerInt>(a, n, out); return;
            case MASK_K11_T18_OPTIMAL:
                extract_kmer_ct_batch<EXT_K11_T18_OPT, KmerInt>(a, n, out); return;
            case MASK_K11_T21_CODING:
                extract_kmer_ct_batch<EXT_K11_T21_COD, KmerInt>(a, n, out); return;
            case MASK_K11_T21_OPTIMAL:
                extract_kmer_ct_batch<EXT_K11_T21_OPT, KmerInt>(a, n, out); return;
            case MASK_K12_T16_CODING:
                extract_kmer_ct_batch<EXT_K12_T16_COD, KmerInt>(a, n, out); return;
            case MASK_K12_T16_OPTIMAL:
                extract_kmer_ct_batch<EXT_K12_T16_OPT, KmerInt>(a, n, out); return;
            case MASK_K12_T18_CODING:
                extract_kmer_ct_batch<EXT_K12_T18_COD, KmerInt>(a, n, out); return;
            case MASK_K12_T18_OPTIMAL:
                extract_kmer_ct_batch<EXT_K12_T18_OPT, KmerInt>(a, n, out); return;
            case MASK_K12_T21_CODING:
                extract_kmer_ct_batch<EXT_K12_T21_COD, KmerInt>(a, n, out); return;
            case MASK_K12_T21_OPTIMAL:
                extract_kmer_ct_batch<EXT_K12_T21_OPT, KmerInt>(a, n, out); return;
            default: fallback(); return;
        }
    } else {
        fallback();
    }
}

template void extract_for_mask_batch<std::uint16_t>(
    const std::uint64_t*, std::size_t, std::uint32_t, std::uint16_t*) noexcept;
template void extract_for_mask_batch<std::uint32_t>(
    const std::uint64_t*, std::size_t, std::uint32_t, std::uint32_t*) noexcept;

} // namespace ikafssn
