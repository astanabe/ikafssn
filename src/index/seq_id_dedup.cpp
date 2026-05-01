// Phase 5i — runtime dispatcher for the per-tier SIMD dedup kernel.
//
// seq_id_dedup_tier.cpp is compiled once per ISA tier as
// ikafssn_dedup_<tier> OBJECT libraries (declared in the top-level
// CMakeLists.txt) under per-tier inner namespaces.  This TU itself is
// built with the project's default flags so it must not include any
// SIMD intrinsics directly — it only routes through function pointers.
//
// Active tiers per architecture are identical to the .kix/.kpx codec
// dispatcher (src/index/pfd_codec.cpp):
//   x86_64 : sse42 / avx2 / avx512bw / avx512vbmi2
//   aarch64: neon
//
// init_simd_dispatch is invoked here on first use (idempotent) so unit
// tests that forget to bootstrap SIMD detection still get a valid
// vtable.

#include "index/seq_id_dedup.hpp"
#include "util/simd_dispatch.hpp"

#include <cstdio>
#include <cstdlib>

#if defined(__x86_64__) || defined(__i386__)
  #define IKAFSSN_DEDUP_HAS_X86  1
#else
  #define IKAFSSN_DEDUP_HAS_X86  0
#endif

#if defined(__aarch64__) || defined(_M_ARM64)
  #define IKAFSSN_DEDUP_HAS_NEON 1
#else
  #define IKAFSSN_DEDUP_HAS_NEON 0
#endif

namespace ikafssn::seq_id_dedup {

#define DECLARE_DEDUP_TIER_NS(ns)                                       \
    namespace ns {                                                      \
        std::uint32_t dedup_seq_ids(const std::uint32_t*, std::uint32_t,\
                                    std::uint32_t*, std::uint32_t*);     \
    }

#if IKAFSSN_DEDUP_HAS_X86
DECLARE_DEDUP_TIER_NS(ikafssn_dedup_sse42)
DECLARE_DEDUP_TIER_NS(ikafssn_dedup_avx2)
DECLARE_DEDUP_TIER_NS(ikafssn_dedup_avx512bw)
DECLARE_DEDUP_TIER_NS(ikafssn_dedup_avx512vbmi2)
#endif

#if IKAFSSN_DEDUP_HAS_NEON
DECLARE_DEDUP_TIER_NS(ikafssn_dedup_neon)
#endif

#undef DECLARE_DEDUP_TIER_NS

namespace {

struct VTable {
    std::uint32_t (*dedup)(const std::uint32_t*, std::uint32_t,
                           std::uint32_t*, std::uint32_t*);
    const char* tier_name;
};

const VTable& active_vtable() {
    static const VTable instance = []() -> VTable {
        init_simd_dispatch(nullptr);
        const SimdCap cap = current_simd_cap();
        (void)cap;

#if IKAFSSN_DEDUP_HAS_X86
        if (cap >= SimdCap::AVX512VBMI2) {
            return { ikafssn_dedup_avx512vbmi2::dedup_seq_ids, "avx512vbmi2" };
        }
        if (cap >= SimdCap::AVX512BW) {
            return { ikafssn_dedup_avx512bw::dedup_seq_ids, "avx512bw" };
        }
        if (cap >= SimdCap::AVX2) {
            return { ikafssn_dedup_avx2::dedup_seq_ids, "avx2" };
        }
        if (cap >= SimdCap::SSE42) {
            return { ikafssn_dedup_sse42::dedup_seq_ids, "sse42" };
        }
        std::fprintf(stderr,
            "ikafssn: seq_id_dedup requires SSE4.2; current CPU tier is below SSE4.2.\n");
        std::exit(2);
#endif // IKAFSSN_DEDUP_HAS_X86

#if IKAFSSN_DEDUP_HAS_NEON
        if (cap >= SimdCap::NEON) {
            return { ikafssn_dedup_neon::dedup_seq_ids, "neon" };
        }
        std::fprintf(stderr,
            "ikafssn: seq_id_dedup requires NEON; current CPU tier is below NEON.\n");
        std::exit(2);
#endif // IKAFSSN_DEDUP_HAS_NEON

#if !IKAFSSN_DEDUP_HAS_X86 && !IKAFSSN_DEDUP_HAS_NEON
        std::fprintf(stderr,
            "ikafssn: seq_id_dedup is not implemented for this architecture.\n");
        std::exit(2);
#endif
    }();
    return instance;
}

} // anonymous namespace

std::uint32_t dedup_seq_ids(const std::uint32_t* sid, std::uint32_t n,
                            std::uint32_t* distinct_sid_out,
                            std::uint32_t* occ_count_out) {
    return active_vtable().dedup(sid, n, distinct_sid_out, occ_count_out);
}

const char* active_tier_name() { return active_vtable().tier_name; }

} // namespace ikafssn::seq_id_dedup
