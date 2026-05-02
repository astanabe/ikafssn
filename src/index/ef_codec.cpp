// Phase 7a — Elias-Fano dictionary codec runtime dispatcher.
//
// In 7a only a single tier exists ("scalar"); all calls route to
// ikafssn::ef::ikafssn_ef_scalar.  Phase 7b will add per-ISA tier
// objects (sse42 / avx2 / avx512bw / avx512vbmi2 / neon) and turn
// this dispatcher into the same SimdCap-cascade pattern as
// pfd_codec.cpp.

#include "index/ef_codec.hpp"
#include "index/ef_format.hpp"

#include <cstdint>
#include <cstring>

#if defined(__unix__) || defined(__APPLE__)
  #include <sys/mman.h>
#endif

namespace ikafssn::ef {

// === Forward declarations of the scalar tier (Phase 7a single tier) ===
namespace ikafssn_ef_scalar {

std::size_t encode_dictionary_ef(const std::uint64_t* offsets,
                                 std::size_t D,
                                 std::uint64_t U_raw,
                                 std::vector<std::uint8_t>& out);

std::uint64_t access_dictionary_ef(const EFHeader& hdr,
                                   const std::uint64_t* lower,
                                   const std::uint64_t* upper,
                                   const std::uint64_t* select,
                                   std::uint32_t i) noexcept;

} // namespace ikafssn_ef_scalar

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
    // Phase 7a: single scalar tier; no SimdCap cascade yet.  7b
    // expands this to mirror pfd_codec.cpp.
    static const VTable instance = []() -> VTable {
        return VTable{
            ikafssn_ef_scalar::encode_dictionary_ef,
            ikafssn_ef_scalar::access_dictionary_ef,
            "scalar",
        };
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
