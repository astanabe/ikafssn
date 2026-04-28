// Phase 5d/5f — runtime dispatcher for the per-tier FastPFor wrappers.
//
// pfd_codec_tier.cpp is compiled four times (once per ISA tier) under the
// ikafssn::pfd::ikafssn_pfd_{sse42,avx2,avx512bw,avx512vbmi2} namespaces.
// Here we expose forward declarations of the per-tier API and pick one
// at first use based on the runtime CPU detection done by simd_dispatch.
//
// This TU itself is built with the project's default flags (no -m...
// pinning), so it must not include any FastPFor header or use any
// SIMD intrinsics directly — it only routes through function pointers.
//
// Phase 5h-prelim: on non-x86_64 builds (aarch64 etc.) the four x86 tier
// OBJECT libraries are not produced (FastPFor's simdpackwithoutmask /
// simdunpack are SSE-only).  The corresponding tier-namespace forward
// declarations and dispatch branches are guarded by `#if defined(__x86_64__)`
// here.  On other architectures `active_vtable()` aborts on first use with
// a clear error; Phase 5h proper will introduce a NEON / SVE2 tier family
// (`ikafssn_pfd_neon` / `ikafssn_pfd_sve2`) that takes those branches.

#include "index/pfd_codec.hpp"
#include "util/simd_dispatch.hpp"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <vector>

#if defined(__x86_64__) || defined(__i386__)
  #define IKAFSSN_PFD_HAS_X86 1
#else
  #define IKAFSSN_PFD_HAS_X86 0
#endif

namespace ikafssn::pfd {

#if IKAFSSN_PFD_HAS_X86
// Forward declarations of every per-tier API.
#define DECLARE_TIER_NS(ns)                                                   \
    namespace ns {                                                            \
        std::size_t encode_posting_kix(const std::uint32_t*, std::uint32_t,   \
                                       std::vector<std::uint8_t>&);           \
        std::size_t encode_posting_kpx(const std::uint32_t*,                  \
                                       const std::uint32_t*,                  \
                                       std::uint32_t, std::uint32_t,          \
                                       std::vector<std::uint8_t>&);           \
        bool open_stream_kix(const std::uint8_t*, std::size_t, StreamCtx&);   \
        bool open_stream_kpx(const std::uint8_t*, std::size_t,                \
                             const std::uint32_t*, std::size_t,               \
                             StreamCtx&);                                     \
    }

DECLARE_TIER_NS(ikafssn_pfd_sse42)
DECLARE_TIER_NS(ikafssn_pfd_avx2)
DECLARE_TIER_NS(ikafssn_pfd_avx512bw)
DECLARE_TIER_NS(ikafssn_pfd_avx512vbmi2)

#undef DECLARE_TIER_NS
#endif // IKAFSSN_PFD_HAS_X86

namespace {

struct VTable {
    std::size_t (*encode_kix)(const std::uint32_t*, std::uint32_t,
                              std::vector<std::uint8_t>&);
    std::size_t (*encode_kpx)(const std::uint32_t*, const std::uint32_t*,
                              std::uint32_t, std::uint32_t,
                              std::vector<std::uint8_t>&);
    bool (*open_kix)(const std::uint8_t*, std::size_t, StreamCtx&);
    bool (*open_kpx)(const std::uint8_t*, std::size_t,
                     const std::uint32_t*, std::size_t,
                     StreamCtx&);
    const char* tier_name;
};

// Resolve the active tier exactly once.  C++11+ guarantees the first
// invocation initialises the local static under a thread-safe one-shot;
// subsequent calls are a single relaxed load.
const VTable& active_vtable() {
    static const VTable instance = []() -> VTable {
        // init_simd_dispatch is idempotent (std::call_once internally) so
        // it is safe to invoke here even when main() has already done so.
        // This makes the codec usable from test binaries that forget to
        // initialise SIMD dispatch explicitly.
        init_simd_dispatch(nullptr);
        const SimdCap cap = current_simd_cap();
        (void)cap;

#if IKAFSSN_PFD_HAS_X86
        // AVX-512 VBMI2 implies VBMI implies BW implies F.  We merge VBMI
        // and VBMI2 into a single top tier here because FastPFor's source
        // contains no byte-level shuffle/permute that VBMI/VBMI2 could
        // exploit; the only difference between the two is compiler auto-
        // vectorization quirks.  Phase 5f removed the standalone VBMI tier
        // from SimdCap, so the highest tier we encounter is AVX512VBMI2.
        if (cap >= SimdCap::AVX512VBMI2) {
            return {
                ikafssn_pfd_avx512vbmi2::encode_posting_kix,
                ikafssn_pfd_avx512vbmi2::encode_posting_kpx,
                ikafssn_pfd_avx512vbmi2::open_stream_kix,
                ikafssn_pfd_avx512vbmi2::open_stream_kpx,
                "avx512vbmi2",
            };
        }
        if (cap >= SimdCap::AVX512BW) {
            return {
                ikafssn_pfd_avx512bw::encode_posting_kix,
                ikafssn_pfd_avx512bw::encode_posting_kpx,
                ikafssn_pfd_avx512bw::open_stream_kix,
                ikafssn_pfd_avx512bw::open_stream_kpx,
                "avx512bw",
            };
        }
        if (cap >= SimdCap::AVX2) {
            return {
                ikafssn_pfd_avx2::encode_posting_kix,
                ikafssn_pfd_avx2::encode_posting_kpx,
                ikafssn_pfd_avx2::open_stream_kix,
                ikafssn_pfd_avx2::open_stream_kpx,
                "avx2",
            };
        }
        if (cap >= SimdCap::SSE42) {
            return {
                ikafssn_pfd_sse42::encode_posting_kix,
                ikafssn_pfd_sse42::encode_posting_kpx,
                ikafssn_pfd_sse42::open_stream_kix,
                ikafssn_pfd_sse42::open_stream_kpx,
                "sse42",
            };
        }

        // SSE4.2 is the x86_64 floor (Phase 5f).  init_simd_dispatch()
        // already rejects pre-SSE4.2 CPUs at startup with exit(2), so the
        // only way to reach this branch is a programmer bug — for example
        // calling the codec on a build_disabled (IKAFSSN_ENABLE_SIMD=0)
        // configuration that bypassed the normal init path.  Abort hard.
        std::fprintf(stderr,
            "ikafssn: pfd codec requires SSE4.2; current CPU tier is below SSE4.2.\n"
            "         (Phase 5f treats SSE4.2 as the x86_64 baseline.)\n");
        std::exit(2);
#else
        // Phase 5h-prelim: non-x86_64 builds (aarch64 etc.) have no
        // FastPFor tier wired up yet.  Tests that do not touch the codec
        // (most of the SIMD-kernel suite) link and run fine because this
        // lambda is only entered on the first .kix/.kpx encode/open call.
        // Real index build / search hits this and aborts with a clear
        // message until Phase 5h proper introduces NEON / SVE2 tiers.
        std::fprintf(stderr,
            "ikafssn: pfd codec is not yet implemented for non-x86_64 architectures.\n"
            "         (Phase 5h-prelim shipped scaffolding; Phase 5h proper will add\n"
            "          NEON / SVE2 PForDelta tier objects.  The SIMD-kernel test suite\n"
            "          still runs natively on aarch64.)\n");
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

std::size_t encode_posting_kpx(const std::uint32_t* sid_array,
                               const std::uint32_t* abs_pos_array,
                               std::uint32_t count,
                               std::uint32_t freq_threshold_part,
                               std::vector<std::uint8_t>& out) {
    return active_vtable().encode_kpx(sid_array, abs_pos_array, count,
                                       freq_threshold_part, out);
}

bool open_stream_kix(const std::uint8_t* posting, std::size_t bytes,
                     StreamCtx& ctx) {
    return active_vtable().open_kix(posting, bytes, ctx);
}

bool open_stream_kpx(const std::uint8_t* posting, std::size_t bytes,
                     const std::uint32_t* sid_stream, std::size_t n_sids,
                     StreamCtx& ctx) {
    return active_vtable().open_kpx(posting, bytes, sid_stream, n_sids, ctx);
}

} // namespace ikafssn::pfd
