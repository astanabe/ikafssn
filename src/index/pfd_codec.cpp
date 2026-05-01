// Phase 5d/5f/5h/5i — runtime dispatcher for the per-tier FastPFor wrappers.
//
// pfd_codec_tier.cpp is compiled once per ISA tier under the
// ikafssn::pfd::ikafssn_pfd_<tier> namespaces.  Here we expose forward
// declarations of the per-tier API and pick one at first use based on the
// runtime CPU detection done by simd_dispatch.
//
// This TU itself is built with the project's default flags (no -m...
// pinning), so it must not include any FastPFor header or use any
// SIMD intrinsics directly — it only routes through function pointers.
//
// Active tiers per architecture:
//   x86_64 : sse42 / avx2 / avx512bw / avx512vbmi2  (Phase 5f 4-tier ladder)
//   aarch64: neon                                   (Phase 5h, single tier)
//
// AArch64 SVE / SVE2 capable CPUs are *not* given a separate tier object:
// the ikafssn_pfd_neon library uses SIMDe (https://github.com/simd-everywhere/simde)
// to translate FastPFor's SSE intrinsics into NEON, and a SVE-native
// bitpacker would be orthogonal hand-coded work.  The dispatcher therefore
// routes any aarch64 SimdCap >= NEON to the single neon tier.

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
            PosDecodeScratch&,                                                \
            std::vector<std::vector<std::uint32_t>>&);                        \
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
                     PosDecodeScratch&,
                     std::vector<std::vector<std::uint32_t>>&);
    const char* tier_name;
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
                "avx512vbmi2",
            };
        }
        if (cap >= SimdCap::AVX512BW) {
            return {
                ikafssn_pfd_avx512bw::encode_posting_kix,
                ikafssn_pfd_avx512bw::encode_posting_kpx,
                ikafssn_pfd_avx512bw::open_stream_kix,
                ikafssn_pfd_avx512bw::open_stream_kpx_for_candidates,
                "avx512bw",
            };
        }
        if (cap >= SimdCap::AVX2) {
            return {
                ikafssn_pfd_avx2::encode_posting_kix,
                ikafssn_pfd_avx2::encode_posting_kpx,
                ikafssn_pfd_avx2::open_stream_kix,
                ikafssn_pfd_avx2::open_stream_kpx_for_candidates,
                "avx2",
            };
        }
        if (cap >= SimdCap::SSE42) {
            return {
                ikafssn_pfd_sse42::encode_posting_kix,
                ikafssn_pfd_sse42::encode_posting_kpx,
                ikafssn_pfd_sse42::open_stream_kix,
                ikafssn_pfd_sse42::open_stream_kpx_for_candidates,
                "sse42",
            };
        }
        std::fprintf(stderr,
            "ikafssn: pfd codec requires SSE4.2; current CPU tier is below SSE4.2.\n"
            "         (Phase 5f treats SSE4.2 as the x86_64 baseline.)\n");
        std::exit(2);
#endif // IKAFSSN_PFD_HAS_X86

#if IKAFSSN_PFD_HAS_NEON
        if (cap >= SimdCap::NEON) {
            return {
                ikafssn_pfd_neon::encode_posting_kix,
                ikafssn_pfd_neon::encode_posting_kpx,
                ikafssn_pfd_neon::open_stream_kix,
                ikafssn_pfd_neon::open_stream_kpx_for_candidates,
                "neon",
            };
        }
        std::fprintf(stderr,
            "ikafssn: pfd codec requires NEON; current CPU tier is below NEON.\n"
            "         (Phase 5h treats NEON as the aarch64 baseline.)\n");
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

const char* active_tier_name() {
    return active_vtable().tier_name;
}

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

bool open_stream_kix(const std::uint8_t* posting, std::size_t bytes,
                     StreamCtx& ctx) {
    return active_vtable().open_kix(posting, bytes, ctx);
}

bool open_stream_kpx_for_candidates(const std::uint8_t* posting, std::size_t bytes,
                                    const std::uint32_t* kix_decoded, std::size_t kix_count,
                                    const std::uint32_t* candidates, std::size_t n_candidates,
                                    PosDecodeScratch& scratch,
                                    std::vector<std::vector<std::uint32_t>>& out) {
    return active_vtable().open_kpx(posting, bytes,
                                    kix_decoded, kix_count,
                                    candidates, n_candidates,
                                    scratch, out);
}

} // namespace ikafssn::pfd
