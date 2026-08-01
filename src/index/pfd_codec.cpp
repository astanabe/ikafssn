// Runtime dispatcher for the per-tier FastPFor wrappers.
//
// pfd_codec_tier.cpp is compiled once per ISA tier under the
// ikafssn::pfd::ikafssn_pfd_<tier> namespaces.  This TU is built with
// the project's default flags (no -m... pinning), so it must not include
// any FastPFor header or use any SIMD intrinsics directly — it only
// routes through function pointers.
//
// Active tiers per architecture:
//   x86_64 : sse42 / avx2 / avx512bw / avx512vbmi2
//   aarch64: neon
//
// On aarch64 the ikafssn_pfd_neon library maps the SSE intrinsics onto
// native NEON through FastPFor's own headers/fastpfor_neon.h; SVE / SVE2
// capable CPUs route to the same neon tier.

#include "index/pfd_codec.hpp"
#include "util/simd_dispatch.hpp"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <vector>

#if defined(__x86_64__) || defined(__i386__)
  #define IKAFSSN_PFD_HAS_X86  1
#else
  #define IKAFSSN_PFD_HAS_X86  0
#endif

#if defined(__aarch64__) || defined(_M_ARM64)
  #define IKAFSSN_PFD_HAS_NEON 1
#else
  #define IKAFSSN_PFD_HAS_NEON 0
#endif

namespace ikafssn::pfd {

// Forward declarations of every per-tier API.  The macro body is unchanged
// across architectures; only the set of namespaces declared differs.
#define DECLARE_TIER_NS(ns)                                                   \
    namespace ns {                                                            \
        std::size_t encode_posting_kix(const std::uint32_t*, std::uint32_t,   \
                                       std::vector<std::uint8_t>&);           \
        std::size_t encode_posting_kpx(const std::uint32_t*,                  \
                                       const std::uint32_t*,                  \
                                       std::uint32_t,                         \
                                       const std::uint32_t*,                  \
                                       std::uint32_t,                         \
                                       std::uint32_t,                         \
                                       std::vector<std::uint8_t>&);           \
        bool open_stream_kix(const std::uint8_t*, std::size_t, StreamCtx&);   \
        bool open_stream_kpx_for_candidates(                                  \
            const std::uint8_t*, std::size_t,                                 \
            const std::uint32_t*, std::size_t,                                \
            const std::uint32_t*, std::size_t,                                \
            PosDecodeScratch&);                                               \
        void popcount_kinds(const std::uint8_t*,                              \
                            std::uint32_t,                                    \
                            std::uint32_t*,                                   \
                            std::uint32_t*,                                   \
                            std::uint32_t*) noexcept;                         \
    }

#if IKAFSSN_PFD_HAS_X86
DECLARE_TIER_NS(ikafssn_pfd_sse42)
DECLARE_TIER_NS(ikafssn_pfd_avx2)
DECLARE_TIER_NS(ikafssn_pfd_avx512bw)
DECLARE_TIER_NS(ikafssn_pfd_avx512vbmi2)
#endif

#if IKAFSSN_PFD_HAS_NEON
DECLARE_TIER_NS(ikafssn_pfd_neon)
#endif

#undef DECLARE_TIER_NS

namespace {

struct VTable {
    std::size_t (*encode_kix)(const std::uint32_t*, std::uint32_t,
                              std::vector<std::uint8_t>&);
    std::size_t (*encode_kpx)(const std::uint32_t*, const std::uint32_t*,
                              std::uint32_t,
                              const std::uint32_t*, std::uint32_t,
                              std::uint32_t,
                              std::vector<std::uint8_t>&);
    bool (*open_kix)(const std::uint8_t*, std::size_t, StreamCtx&);
    bool (*open_kpx)(const std::uint8_t*, std::size_t,
                     const std::uint32_t*, std::size_t,
                     const std::uint32_t*, std::size_t,
                     PosDecodeScratch&);
    void (*popcount)(const std::uint8_t*, std::uint32_t,
                     std::uint32_t*, std::uint32_t*, std::uint32_t*) noexcept;
};

const VTable& active_vtable() {
    static const VTable instance = []() -> VTable {
        init_simd_dispatch(nullptr);
        const SimdCap cap = current_simd_cap();
        (void)cap;

#if IKAFSSN_PFD_HAS_X86
        if (cap >= SimdCap::AVX512VBMI2) {
            return {
                ikafssn_pfd_avx512vbmi2::encode_posting_kix,
                ikafssn_pfd_avx512vbmi2::encode_posting_kpx,
                ikafssn_pfd_avx512vbmi2::open_stream_kix,
                ikafssn_pfd_avx512vbmi2::open_stream_kpx_for_candidates,
                ikafssn_pfd_avx512vbmi2::popcount_kinds,
            };
        }
        if (cap >= SimdCap::AVX512BW) {
            return {
                ikafssn_pfd_avx512bw::encode_posting_kix,
                ikafssn_pfd_avx512bw::encode_posting_kpx,
                ikafssn_pfd_avx512bw::open_stream_kix,
                ikafssn_pfd_avx512bw::open_stream_kpx_for_candidates,
                ikafssn_pfd_avx512bw::popcount_kinds,
            };
        }
        if (cap >= SimdCap::AVX2) {
            return {
                ikafssn_pfd_avx2::encode_posting_kix,
                ikafssn_pfd_avx2::encode_posting_kpx,
                ikafssn_pfd_avx2::open_stream_kix,
                ikafssn_pfd_avx2::open_stream_kpx_for_candidates,
                ikafssn_pfd_avx2::popcount_kinds,
            };
        }
        if (cap >= SimdCap::SSE42) {
            return {
                ikafssn_pfd_sse42::encode_posting_kix,
                ikafssn_pfd_sse42::encode_posting_kpx,
                ikafssn_pfd_sse42::open_stream_kix,
                ikafssn_pfd_sse42::open_stream_kpx_for_candidates,
                ikafssn_pfd_sse42::popcount_kinds,
            };
        }
        std::fprintf(stderr,
            "ikafssn: pfd codec requires SSE4.2 (the x86_64 baseline).\n");
        std::exit(2);
#endif // IKAFSSN_PFD_HAS_X86

#if IKAFSSN_PFD_HAS_NEON
        if (cap >= SimdCap::NEON) {
            return {
                ikafssn_pfd_neon::encode_posting_kix,
                ikafssn_pfd_neon::encode_posting_kpx,
                ikafssn_pfd_neon::open_stream_kix,
                ikafssn_pfd_neon::open_stream_kpx_for_candidates,
                ikafssn_pfd_neon::popcount_kinds,
            };
        }
        std::fprintf(stderr,
            "ikafssn: pfd codec requires NEON (the aarch64 baseline).\n");
        std::exit(2);
#endif // IKAFSSN_PFD_HAS_NEON

#if !IKAFSSN_PFD_HAS_X86 && !IKAFSSN_PFD_HAS_NEON
        std::fprintf(stderr,
            "ikafssn: pfd codec is not implemented for this architecture.\n"
            "         (Only x86_64 and aarch64 are supported.)\n");
        std::exit(2);
#endif
    }();
    return instance;
}

} // anonymous namespace

std::size_t encode_posting_kix(const std::uint32_t* delta_array,
                               std::uint32_t count,
                               std::vector<std::uint8_t>& out) {
    return active_vtable().encode_kix(delta_array, count, out);
}

std::size_t encode_posting_kpx(const std::uint32_t* distinct_sid,
                               const std::uint32_t* occ_count,
                               std::uint32_t distinct_count,
                               const std::uint32_t* abs_pos_array,
                               std::uint32_t position_count,
                               std::uint32_t freq_threshold_part,
                               std::vector<std::uint8_t>& out) {
    return active_vtable().encode_kpx(distinct_sid, occ_count, distinct_count,
                                       abs_pos_array, position_count,
                                       freq_threshold_part, out);
}

bool open_stream_kix(const std::uint8_t* posting_list, std::size_t bytes,
                     StreamCtx& ctx) {
    return active_vtable().open_kix(posting_list, bytes, ctx);
}

bool open_stream_kpx_for_candidates(const std::uint8_t* posting_list, std::size_t bytes,
                                    const std::uint32_t* kix_decoded, std::size_t kix_count,
                                    const std::uint32_t* candidates, std::size_t n_candidates,
                                    PosDecodeScratch& scratch) {
    return active_vtable().open_kpx(posting_list, bytes,
                                    kix_decoded, kix_count,
                                    candidates, n_candidates,
                                    scratch);
}

void popcount_kinds(const std::uint8_t* km, std::uint32_t distinct_count,
                    std::uint32_t* p_partition, std::uint32_t* p_short1,
                    std::uint32_t* p_short2) noexcept {
    active_vtable().popcount(km, distinct_count, p_partition, p_short1, p_short2);
}

} // namespace ikafssn::pfd
