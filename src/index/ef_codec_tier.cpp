// Phase 7a — Elias-Fano codec, scalar reference body.
//
// Phase 7b will add per-ISA tier specialisations under
// ikafssn_ef_{sse42,avx2,avx512bw,avx512vbmi2,neon} following the
// FastPFor / dedup ladder pattern.  For 7a this TU compiles once
// under default flags (no -m... pinning) and exposes the body in
// the ikafssn::ef::ikafssn_ef_scalar namespace, which the dispatcher
// in ef_codec.cpp routes to unconditionally.

#include "index/ef_codec.hpp"
#include "index/ef_format.hpp"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <vector>

#if defined(__unix__) || defined(__APPLE__)
  #include <sys/mman.h>
#endif

namespace ikafssn::ef::ikafssn_ef_scalar {

namespace {

// floor(log2(n)) for n >= 1.  Returns 0 when n == 0.
inline std::uint8_t floor_log2_u64(std::uint64_t n) noexcept {
    if (n == 0) return 0;
    return static_cast<std::uint8_t>(63 - __builtin_clzll(n));
}

inline std::uint64_t mask_for_bits(std::uint8_t l) noexcept {
    return (l >= 64) ? ~std::uint64_t{0}
                     : ((std::uint64_t{1} << l) - std::uint64_t{1});
}

// Write `l` bits of `value` into `lower` starting at bit position
// `bit_pos`.  Caller guarantees `lower` has room for one possible
// spill word past the high end.
inline void write_lower_bits(std::uint64_t* lower,
                             std::uint64_t bit_pos,
                             std::uint64_t value,
                             std::uint8_t l) noexcept {
    if (l == 0) return;
    std::uint64_t word_idx     = bit_pos >> 6;
    std::uint64_t bit_in_word  = bit_pos & 63;
    lower[word_idx] |= value << bit_in_word;
    std::uint64_t spill = bit_in_word + l;
    if (spill > 64) {
        lower[word_idx + 1] |= value >> (64 - bit_in_word);
    }
}

inline std::uint64_t read_lower_bits(const std::uint64_t* lower,
                                     std::uint64_t bit_pos,
                                     std::uint8_t l,
                                     std::uint64_t mask_l) noexcept {
    if (l == 0) return 0;
    std::uint64_t word_idx    = bit_pos >> 6;
    std::uint64_t bit_in_word = bit_pos & 63;
    std::uint64_t lo = (lower[word_idx] >> bit_in_word) & mask_l;
    if (bit_in_word + l > 64) {
        lo |= (lower[word_idx + 1] << (64 - bit_in_word)) & mask_l;
    }
    return lo;
}

// Locate the bit position of the j-th (0-indexed) set bit AT OR AFTER
// upper-bit position `start_bit_pos`, with `j == 0` returning the
// first set bit at or after `start_bit_pos`.  Scalar word-by-word
// scan with __builtin_popcountll + __builtin_ctzll.
inline std::uint64_t select1_after(const std::uint64_t* upper,
                                   std::uint64_t upper_bits_total,
                                   std::uint64_t start_bit_pos,
                                   std::uint64_t j) noexcept {
    std::uint64_t pos = start_bit_pos;
    std::uint64_t found = 0;
    while (pos < upper_bits_total) {
        std::uint64_t word_idx    = pos >> 6;
        std::uint64_t bit_in_word = pos & 63;
        std::uint64_t w = upper[word_idx] >> bit_in_word;
        std::uint64_t cnt = static_cast<std::uint64_t>(__builtin_popcountll(w));
        if (found + cnt > j) {
            std::uint64_t need = j - found;  // 0-indexed bit within w
            // Pop low bits until we reach the (need)-th set bit.
            for (std::uint64_t k = 0; k < need; ++k) {
                w &= w - 1;
            }
            return pos + static_cast<std::uint64_t>(__builtin_ctzll(w));
        }
        found += cnt;
        pos = (word_idx + 1) << 6;  // advance to next word boundary
    }
    return upper_bits_total;  // not found (caller error)
}

} // anonymous namespace

std::size_t encode_dictionary_ef(const std::uint64_t* offsets,
                                 std::size_t D,
                                 std::uint64_t U_raw,
                                 std::vector<std::uint8_t>& out) {
    std::size_t out_start = out.size();

    EFHeader hdr{};
    std::memcpy(hdr.magic, EF_MAGIC, 4);
    hdr.D = static_cast<std::uint64_t>(D);

    if (D == 0) {
        hdr.l = 0;
        hdr.U = 0;
        hdr.upper_bits = 0;
        hdr.select_count = 0;
        out.resize(out_start + sizeof(EFHeader));
        std::memcpy(out.data() + out_start, &hdr, sizeof(EFHeader));
        return sizeof(EFHeader);
    }

    // Strict-monotonic universe: ef[i] = offsets[i] + i ∈ [0, U_raw + D)
    std::uint64_t U_ef = U_raw + static_cast<std::uint64_t>(D);
    std::uint8_t l = (U_ef > static_cast<std::uint64_t>(D))
                         ? floor_log2_u64(U_ef / D)
                         : std::uint8_t{0};
    if (l > 63) l = 63;
    std::uint64_t mask_l = mask_for_bits(l);

    std::uint64_t max_ef = offsets[D - 1] + (D - 1);
    std::uint64_t max_high = max_ef >> l;
    std::uint64_t upper_bits_total = static_cast<std::uint64_t>(D) + max_high + 1;

    std::uint64_t lower_bits_total = static_cast<std::uint64_t>(D) * l;
    // Pad to whole u64 words; write_lower_bits's spill needs one extra.
    std::size_t lower_words = static_cast<std::size_t>((lower_bits_total + 63) / 64);
    std::size_t upper_words = static_cast<std::size_t>((upper_bits_total + 63) / 64);

    std::uint32_t select_count = static_cast<std::uint32_t>(
        (static_cast<std::uint64_t>(D) + kSelectStep - 1) / kSelectStep);

    hdr.l = l;
    hdr.U = U_ef;
    hdr.upper_bits = static_cast<std::uint32_t>(upper_bits_total);
    hdr.select_count = select_count;

    std::size_t blob_bytes = sizeof(EFHeader)
                           + lower_words * 8
                           + upper_words * 8
                           + static_cast<std::size_t>(select_count) * 8;
    out.resize(out_start + blob_bytes);
    std::uint8_t* base = out.data() + out_start;
    std::memcpy(base, &hdr, sizeof(EFHeader));

    std::uint64_t* lower = reinterpret_cast<std::uint64_t*>(base + sizeof(EFHeader));
    std::uint64_t* upper = lower + lower_words;
    std::uint64_t* sel   = upper + upper_words;
    std::memset(lower, 0, lower_words * 8);
    std::memset(upper, 0, upper_words * 8);

    for (std::size_t i = 0; i < D; ++i) {
        std::uint64_t v = offsets[i] + static_cast<std::uint64_t>(i);
        std::uint64_t lo = v & mask_l;
        std::uint64_t hi = v >> l;
        std::uint64_t up_pos = hi + static_cast<std::uint64_t>(i);

        write_lower_bits(lower, static_cast<std::uint64_t>(i) * l, lo, l);
        upper[up_pos >> 6] |= std::uint64_t{1} << (up_pos & 63);

        if (i % kSelectStep == 0) {
            sel[i / kSelectStep] = up_pos;
        }
    }

    return blob_bytes;
}

bool open_dictionary_ef(EFDictionary& self,
                        const std::uint8_t* data,
                        std::size_t bytes,
                        std::size_t& blob_bytes_out,
                        const std::uint64_t*& lower_out,
                        const std::uint64_t*& upper_out,
                        const std::uint64_t*& select_out) noexcept {
    if (bytes < sizeof(EFHeader)) return false;
    EFHeader hdr;
    std::memcpy(&hdr, data, sizeof(EFHeader));
    if (std::memcmp(hdr.magic, EF_MAGIC, 4) != 0) return false;
    if (hdr.l > 63) return false;

    std::size_t lower_words = static_cast<std::size_t>(
        (hdr.D * static_cast<std::uint64_t>(hdr.l) + 63) / 64);
    std::size_t upper_words = static_cast<std::size_t>(
        (static_cast<std::uint64_t>(hdr.upper_bits) + 63) / 64);
    std::size_t blob_bytes = sizeof(EFHeader)
                           + lower_words * 8
                           + upper_words * 8
                           + static_cast<std::size_t>(hdr.select_count) * 8;
    if (bytes < blob_bytes) return false;

    blob_bytes_out = blob_bytes;
    const std::uint8_t* p = data + sizeof(EFHeader);
    lower_out  = reinterpret_cast<const std::uint64_t*>(p);
    p += lower_words * 8;
    upper_out  = reinterpret_cast<const std::uint64_t*>(p);
    p += upper_words * 8;
    select_out = reinterpret_cast<const std::uint64_t*>(p);

    self = EFDictionary{};  // re-init; caller populates remaining members
    (void)hdr;
    return true;
}

std::uint64_t access_dictionary_ef(const EFHeader& hdr,
                                   const std::uint64_t* lower,
                                   const std::uint64_t* upper,
                                   const std::uint64_t* select,
                                   std::uint32_t i) noexcept {
    const std::uint8_t l = hdr.l;
    const std::uint64_t mask_l = mask_for_bits(l);
    const std::uint64_t upper_bits_total = static_cast<std::uint64_t>(hdr.upper_bits);

    std::uint32_t sample_idx = i / kSelectStep;
    std::uint32_t base_rank  = sample_idx * kSelectStep;
    std::uint64_t base_pos   = select[sample_idx];
    std::uint64_t up_pos;
    if (i == base_rank) {
        up_pos = base_pos;
    } else {
        up_pos = select1_after(upper,
                               upper_bits_total,
                               base_pos + 1,
                               static_cast<std::uint64_t>(i - base_rank - 1));
    }

    std::uint64_t hi = up_pos - static_cast<std::uint64_t>(i);
    std::uint64_t lo = read_lower_bits(lower,
                                       static_cast<std::uint64_t>(i) * l,
                                       l, mask_l);
    std::uint64_t ef_value = (hi << l) | lo;
    return ef_value - static_cast<std::uint64_t>(i);
}

} // namespace ikafssn::ef::ikafssn_ef_scalar
