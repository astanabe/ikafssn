// Phase 7b — Elias-Fano dictionary codec runtime dispatcher.
//
// Mirrors the FastPFor / dedup ladder in pfd_codec.cpp:
//   x86_64 : sse42 / avx2 / avx512bw / avx512vbmi2  (SSE4.2 floor)
//   aarch64: neon                                    (NEON floor)
//
// CPUs below SSE4.2 (x86) or NEON (aarch64) are rejected at startup
// with exit(2), matching the rest of the codebase's policy.

#include "index/ef_codec.hpp"
#include "index/ef_format.hpp"
#include "util/simd_dispatch.hpp"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#if defined(__unix__) || defined(__APPLE__)
  #include <sys/mman.h>
#endif

#if defined(__x86_64__) || defined(__i386__)
  #define IKAFSSN_EF_HAS_X86  1
#else
  #define IKAFSSN_EF_HAS_X86  0
#endif

#if defined(__aarch64__) || defined(_M_ARM64)
  #define IKAFSSN_EF_HAS_NEON 1
#else
  #define IKAFSSN_EF_HAS_NEON 0
#endif

namespace ikafssn::ef {

#define DECLARE_EF_TIER_NS(ns)                                                 \
    namespace ns {                                                             \
        std::size_t encode_dictionary_ef(const std::uint64_t*,                 \
                                         std::size_t,                          \
                                         std::uint64_t,                        \
                                         std::vector<std::uint8_t>&);          \
        std::uint64_t access_dictionary_ef(const EFHeader&,                    \
                                           const std::uint64_t*,               \
                                           const std::uint64_t*,               \
                                           const std::uint64_t*,               \
                                           std::uint32_t) noexcept;            \
    }

#if IKAFSSN_EF_HAS_X86
DECLARE_EF_TIER_NS(ikafssn_ef_sse42)
DECLARE_EF_TIER_NS(ikafssn_ef_avx2)
DECLARE_EF_TIER_NS(ikafssn_ef_avx512bw)
DECLARE_EF_TIER_NS(ikafssn_ef_avx512vbmi2)
#endif
#if IKAFSSN_EF_HAS_NEON
DECLARE_EF_TIER_NS(ikafssn_ef_neon)
#endif

#undef DECLARE_EF_TIER_NS

namespace {

struct VTable {
    std::size_t (*encode_ef)(const std::uint64_t*, std::size_t, std::uint64_t,
                             std::vector<std::uint8_t>&);
    std::uint64_t (*access_ef)(const EFHeader&,
                               const std::uint64_t*, const std::uint64_t*,
                               const std::uint64_t*, std::uint32_t) noexcept;
    const char* tier_name;
};

const VTable& active_vtable() {
    static const VTable instance = []() -> VTable {
        init_simd_dispatch(nullptr);
        const SimdCap cap = current_simd_cap();
        (void)cap;

#if IKAFSSN_EF_HAS_X86
        if (cap >= SimdCap::AVX512VBMI2) {
            return {
                ikafssn_ef_avx512vbmi2::encode_dictionary_ef,
                ikafssn_ef_avx512vbmi2::access_dictionary_ef,
                "avx512vbmi2",
            };
        }
        if (cap >= SimdCap::AVX512BW) {
            return {
                ikafssn_ef_avx512bw::encode_dictionary_ef,
                ikafssn_ef_avx512bw::access_dictionary_ef,
                "avx512bw",
            };
        }
        if (cap >= SimdCap::AVX2) {
            return {
                ikafssn_ef_avx2::encode_dictionary_ef,
                ikafssn_ef_avx2::access_dictionary_ef,
                "avx2",
            };
        }
        if (cap >= SimdCap::SSE42) {
            return {
                ikafssn_ef_sse42::encode_dictionary_ef,
                ikafssn_ef_sse42::access_dictionary_ef,
                "sse42",
            };
        }
        std::fprintf(stderr,
            "ikafssn: ef codec requires SSE4.2; current CPU tier is below SSE4.2.\n"
            "         (Phase 5f treats SSE4.2 as the x86_64 baseline.)\n");
        std::exit(2);
#endif

#if IKAFSSN_EF_HAS_NEON
        if (cap >= SimdCap::NEON) {
            return {
                ikafssn_ef_neon::encode_dictionary_ef,
                ikafssn_ef_neon::access_dictionary_ef,
                "neon",
            };
        }
        std::fprintf(stderr,
            "ikafssn: ef codec requires NEON; current CPU tier is below NEON.\n"
            "         (Phase 5h treats NEON as the aarch64 baseline.)\n");
        std::exit(2);
#endif

#if !IKAFSSN_EF_HAS_X86 && !IKAFSSN_EF_HAS_NEON
        std::fprintf(stderr,
            "ikafssn: ef codec is not implemented for this architecture.\n"
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

std::size_t encode_dictionary_ef(const std::uint64_t* offsets,
                                 std::size_t D,
                                 std::uint64_t U_raw,
                                 std::vector<std::uint8_t>& out) {
    return active_vtable().encode_ef(offsets, D, U_raw, out);
}

bool EFDictionary::open(const std::uint8_t* data, std::size_t bytes) {
    if (bytes < sizeof(EFHeader)) return false;
    EFHeader hdr;
    std::memcpy(&hdr, data, sizeof(EFHeader));
    if (std::memcmp(hdr.magic, EF_MAGIC, 4) != 0) return false;
    if (hdr.l > 63) return false;

    const std::uint64_t lower_bits_total =
        hdr.D * static_cast<std::uint64_t>(hdr.l);
    const std::size_t lower_words = static_cast<std::size_t>((lower_bits_total + 63) / 64);
    const std::size_t upper_words = static_cast<std::size_t>(
        (static_cast<std::uint64_t>(hdr.upper_bits) + 63) / 64);
    const std::size_t need = sizeof(EFHeader)
                           + lower_words * 8
                           + upper_words * 8
                           + static_cast<std::size_t>(hdr.select_count) * 8;
    if (bytes < need) return false;

    header_ = hdr;
    data_ = data;
    blob_bytes_ = need;

    const std::uint8_t* p = data + sizeof(EFHeader);
    lower_  = reinterpret_cast<const std::uint64_t*>(p);
    p += lower_words * 8;
    upper_  = reinterpret_cast<const std::uint64_t*>(p);
    p += upper_words * 8;
    select_ = reinterpret_cast<const std::uint64_t*>(p);

    D_ = static_cast<std::size_t>(hdr.D);
    l_ = hdr.l;
    mask_l_ = (l_ >= 64) ? ~std::uint64_t{0}
                         : ((std::uint64_t{1} << l_) - std::uint64_t{1});
    upper_bits_total_ = static_cast<std::uint64_t>(hdr.upper_bits);
    select_count_ = hdr.select_count;
    return true;
}

void EFDictionary::apply_madvise(bool willneed) const noexcept {
#if defined(__unix__) || defined(__APPLE__)
    if (data_ == nullptr || blob_bytes_ == 0) return;
    int advice = willneed ? POSIX_MADV_WILLNEED : POSIX_MADV_RANDOM;
    posix_madvise(const_cast<void*>(static_cast<const void*>(data_)),
                  blob_bytes_, advice);
#else
    (void)willneed;
#endif
}

std::uint64_t EFDictionary::access(std::uint32_t i) const noexcept {
    if (i >= D_) return 0;
    return active_vtable().access_ef(header_, lower_, upper_, select_, i);
}

void EFDictionary::access_pair(std::uint32_t i,
                               std::uint64_t& start,
                               std::uint64_t& end) const noexcept {
    start = access(i);
    end = access(i + 1);
}

} // namespace ikafssn::ef
