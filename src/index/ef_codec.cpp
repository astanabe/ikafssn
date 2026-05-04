// Phase 7b — Elias-Fano dictionary codec runtime dispatcher.
//
// Mirrors the FastPFor / dedup ladder in pfd_codec.cpp:
//   x86_64 : sse42 / avx2 / avx512bw / avx512vbmi2  (SSE4.2 floor)
//   aarch64: neon                                    (NEON floor)
//
// CPUs below SSE4.2 (x86) or NEON (aarch64) are rejected at startup
// with exit(2), matching the rest of the codebase's policy.
//
// Encoding is implemented here (non-tier) — the encoder body is pure
// scalar arithmetic so per-tier compilation buys nothing, and putting
// it in the dispatcher TU lets the chunk-parallel orchestration use
// TBB without dragging it into every per-tier OBJECT library.

#include "index/ef_codec.hpp"
#include "index/ef_format.hpp"
#include "util/simd_dispatch.hpp"

#include <atomic>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/task_arena.h>

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
                ikafssn_ef_avx512vbmi2::access_dictionary_ef,
                "avx512vbmi2",
            };
        }
        if (cap >= SimdCap::AVX512BW) {
            return {
                ikafssn_ef_avx512bw::access_dictionary_ef,
                "avx512bw",
            };
        }
        if (cap >= SimdCap::AVX2) {
            return {
                ikafssn_ef_avx2::access_dictionary_ef,
                "avx2",
            };
        }
        if (cap >= SimdCap::SSE42) {
            return {
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

inline std::uint8_t floor_log2_u64(std::uint64_t n) noexcept {
    if (n == 0) return 0;
    return static_cast<std::uint8_t>(63 - __builtin_clzll(n));
}

inline std::uint64_t mask_for_bits(std::uint8_t l) noexcept {
    return (l >= 64) ? ~std::uint64_t{0}
                     : ((std::uint64_t{1} << l) - std::uint64_t{1});
}

// Below this threshold the parallel orchestration overhead dominates
// the actual encode work; fall back to a single-thread loop.
constexpr std::size_t kEFParallelThreshold = 1u << 14;  // 16384 entries

} // anonymous namespace

const char* active_tier_name() {
    return active_vtable().tier_name;
}

// Encode the EF dictionary blob.  This is non-tier (the body is pure
// scalar arithmetic, so SIMD-level specialisation buys nothing) and
// runs the per-entry loop chunk-parallel via TBB.
//
// Each chunk owns a slice of the input [i_begin, i_end) and writes its
// lower-bit and upper-bit contributions into a chunk-local buffer that
// covers exactly the global word range it touches.  Lower-bit ranges
// are bit-contiguous (chunk c+1 starts at the bit immediately after
// chunk c's last) so adjacent chunks share at most one global word at
// the boundary.  Upper-bit positions are monotonic in i (offsets[i] is
// non-decreasing, so up_pos = ((offsets[i]+i)>>l)+i is non-decreasing
// too), which gives the same one-word boundary-overlap property.
//
// The merge therefore writes inner words with plain stores (each is
// touched by exactly one chunk) and ORs boundary words with
// std::atomic_ref to handle the at-most-one-word overlap with each
// neighbouring chunk safely without serialising the whole merge.
std::size_t encode_dictionary_ef(const std::uint64_t* offsets,
                                 std::size_t D,
                                 std::uint64_t U_raw,
                                 std::vector<std::uint8_t>& out) {
    const std::size_t out_start = out.size();

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

    const std::uint64_t U_ef = U_raw + static_cast<std::uint64_t>(D);
    std::uint8_t l = (U_ef > static_cast<std::uint64_t>(D))
                         ? floor_log2_u64(U_ef / D)
                         : std::uint8_t{0};
    if (l > 63) l = 63;
    const std::uint64_t mask_l = mask_for_bits(l);

    const std::uint64_t max_ef = offsets[D - 1] + (D - 1);
    const std::uint64_t max_high = max_ef >> l;
    const std::uint64_t upper_bits_total = static_cast<std::uint64_t>(D) + max_high + 1;

    const std::uint64_t lower_bits_total = static_cast<std::uint64_t>(D) * l;
    const std::size_t lower_words = static_cast<std::size_t>((lower_bits_total + 63) / 64);
    const std::size_t upper_words = static_cast<std::size_t>((upper_bits_total + 63) / 64);

    const std::uint32_t select_count = static_cast<std::uint32_t>(
        (static_cast<std::uint64_t>(D) + kSelectStep - 1) / kSelectStep);

    hdr.l = l;
    hdr.U = U_ef;
    hdr.upper_bits = static_cast<std::uint32_t>(upper_bits_total);
    hdr.select_count = select_count;

    const std::size_t blob_bytes = sizeof(EFHeader)
                                 + lower_words * 8
                                 + upper_words * 8
                                 + static_cast<std::size_t>(select_count) * 8;
    out.resize(out_start + blob_bytes);
    std::uint8_t* base = out.data() + out_start;
    std::memcpy(base, &hdr, sizeof(EFHeader));

    std::uint64_t* const lower = reinterpret_cast<std::uint64_t*>(base + sizeof(EFHeader));
    std::uint64_t* const upper = lower + lower_words;
    std::uint64_t* const sel   = upper + upper_words;
    std::memset(lower, 0, lower_words * 8);
    std::memset(upper, 0, upper_words * 8);

    // Small inputs: the parallel orchestration overhead exceeds the
    // serial cost.  Run a tight single-thread loop and bail out.
    if (D < kEFParallelThreshold) {
        for (std::size_t i = 0; i < D; ++i) {
            const std::uint64_t v = offsets[i] + static_cast<std::uint64_t>(i);
            const std::uint64_t lo = v & mask_l;
            const std::uint64_t hi = v >> l;
            const std::uint64_t up_pos = hi + static_cast<std::uint64_t>(i);

            if (l > 0) {
                const std::uint64_t bit_pos = static_cast<std::uint64_t>(i) * l;
                const std::uint64_t word_idx = bit_pos >> 6;
                const std::uint64_t bit_in_word = bit_pos & 63;
                lower[word_idx] |= lo << bit_in_word;
                if (bit_in_word + l > 64) {
                    lower[word_idx + 1] |= lo >> (64 - bit_in_word);
                }
            }
            upper[up_pos >> 6] |= std::uint64_t{1} << (up_pos & 63);
            if (i % kSelectStep == 0) {
                sel[i / kSelectStep] = up_pos;
            }
        }
        return blob_bytes;
    }

    // Pick chunk count: cap at concurrency, ensure each chunk has a
    // worthwhile slice (>= ~16k entries).
    std::size_t num_chunks = static_cast<std::size_t>(
        tbb::this_task_arena::max_concurrency());
    if (num_chunks == 0) num_chunks = 1;
    const std::size_t min_chunk_entries = 1u << 14;
    const std::size_t max_chunks_by_size = (D + min_chunk_entries - 1) / min_chunk_entries;
    if (num_chunks > max_chunks_by_size) num_chunks = max_chunks_by_size;
    if (num_chunks == 0) num_chunks = 1;

    std::vector<std::size_t> i_bounds(num_chunks + 1);
    for (std::size_t c = 0; c <= num_chunks; ++c) {
        i_bounds[c] = (D * c) / num_chunks;
    }

    struct ChunkOut {
        std::vector<std::uint64_t> local_lower;
        std::vector<std::uint64_t> local_upper;
        std::size_t lower_word_begin = 0;
        std::size_t upper_word_begin = 0;
    };
    std::vector<ChunkOut> chunks(num_chunks);

    tbb::parallel_for(std::size_t{0}, num_chunks, [&](std::size_t c) {
        const std::size_t ib = i_bounds[c];
        const std::size_t ie = i_bounds[c + 1];
        if (ib >= ie) return;

        // Lower-bit window: bit positions [ib*l, ie*l).
        const std::uint64_t lb_bit_begin = static_cast<std::uint64_t>(ib) * l;
        const std::uint64_t lb_bit_end   = static_cast<std::uint64_t>(ie) * l;
        const std::size_t lb_word_begin = static_cast<std::size_t>(lb_bit_begin >> 6);
        // (lb_bit_end + 63) / 64; if l == 0 this collapses to lb_word_begin
        // with zero words, which is fine (the inner loop skips lower writes).
        const std::size_t lb_word_end =
            (l > 0) ? static_cast<std::size_t>((lb_bit_end + 63) >> 6)
                    : lb_word_begin;

        // Upper-bit window: derive from monotonic up_pos at first / last i.
        const std::uint64_t up_first = ((offsets[ib] + ib) >> l) + ib;
        const std::uint64_t up_last  = ((offsets[ie - 1] + (ie - 1)) >> l) + (ie - 1);
        const std::size_t ub_word_begin = static_cast<std::size_t>(up_first >> 6);
        const std::size_t ub_word_end   = static_cast<std::size_t>((up_last >> 6) + 1);

        ChunkOut& co = chunks[c];
        co.lower_word_begin = lb_word_begin;
        co.upper_word_begin = ub_word_begin;
        co.local_lower.assign(lb_word_end - lb_word_begin, 0);
        co.local_upper.assign(ub_word_end - ub_word_begin, 0);

        std::uint64_t* const lo_buf = co.local_lower.data();
        std::uint64_t* const up_buf = co.local_upper.data();

        for (std::size_t i = ib; i < ie; ++i) {
            const std::uint64_t v = offsets[i] + static_cast<std::uint64_t>(i);
            const std::uint64_t lo = v & mask_l;
            const std::uint64_t hi = v >> l;
            const std::uint64_t up_pos = hi + static_cast<std::uint64_t>(i);

            if (l > 0) {
                const std::uint64_t bit_pos = static_cast<std::uint64_t>(i) * l;
                const std::size_t local_word =
                    static_cast<std::size_t>(bit_pos >> 6) - lb_word_begin;
                const std::uint64_t bit_in_word = bit_pos & 63;
                lo_buf[local_word] |= lo << bit_in_word;
                if (bit_in_word + l > 64) {
                    lo_buf[local_word + 1] |= lo >> (64 - bit_in_word);
                }
            }
            const std::size_t up_local =
                static_cast<std::size_t>(up_pos >> 6) - ub_word_begin;
            up_buf[up_local] |= std::uint64_t{1} << (up_pos & 63);

            // Each i maps to a distinct sample slot, so writes are
            // race-free across chunks without atomics.
            if (i % kSelectStep == 0) {
                sel[i / kSelectStep] = up_pos;
            }
        }
    });

    // Merge.  Inner words are owned by exactly one chunk (plain store).
    // Boundary words may overlap with the previous / next chunk by at
    // most one word, so OR them through atomic_ref.
    tbb::parallel_for(std::size_t{0}, num_chunks, [&](std::size_t c) {
        const ChunkOut& co = chunks[c];
        const std::size_t lo_n = co.local_lower.size();
        for (std::size_t w = 0; w < lo_n; ++w) {
            const std::uint64_t v = co.local_lower[w];
            if (v == 0) continue;
            const std::size_t global_w = co.lower_word_begin + w;
            const bool boundary = (w == 0) || (w + 1 == lo_n);
            if (boundary) {
                std::atomic_ref<std::uint64_t>(lower[global_w])
                    .fetch_or(v, std::memory_order_relaxed);
            } else {
                lower[global_w] |= v;
            }
        }
        const std::size_t up_n = co.local_upper.size();
        for (std::size_t w = 0; w < up_n; ++w) {
            const std::uint64_t v = co.local_upper[w];
            if (v == 0) continue;
            const std::size_t global_w = co.upper_word_begin + w;
            const bool boundary = (w == 0) || (w + 1 == up_n);
            if (boundary) {
                std::atomic_ref<std::uint64_t>(upper[global_w])
                    .fetch_or(v, std::memory_order_relaxed);
            } else {
                upper[global_w] |= v;
            }
        }
    });

    return blob_bytes;
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
