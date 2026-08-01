// Per-tier custom block codec.
//
// This translation unit is compiled four times on x86_64 (once per ISA
// tier) and once on aarch64 via the ikafssn_pfd_<tier> OBJECT libraries
// declared in the top-level CMakeLists.txt.  Each compilation is given:
//
//   -DFastPForLib=FastPForLib_<tier>     renames FastPFor's namespace
//                                         at preprocessor time so the
//                                         per-tier sets of symbols do
//                                         not collide at link time
//   -DIKAFSSN_PFD_TIER_NAME=<tier>       names the tier-specific
//                                         ikafssn::pfd inner namespace
//   -m<arch> ...                          controls the instructions the
//                                         bitpacker primitives and the
//                                         surrounding scalar code in
//                                         this file are allowed to use
//
// For .kix posting lists (sorted distinct seq_id delta stream) we drive
// FastPFor's CompositeCodec<SIMDFastPFor<4>, VariableByte> directly.
//
// For .kpx posting lists (absolute position stream) we keep the per-
// (kmer, seq_id) partition grouping and classify every distinct
// seq_id via a 2-bit kind map.  The short bucket is split into occ=1
// and occ>=2 sub-buckets so single-position clusters drop the u8
// occ_count entirely.  Each FOR-block header is 8 B and stream tails
// use a packed bit-width encoding.

#include "index/pfd_codec.hpp"

// FastPForLib is rewritten by the build system to FastPForLib_<tier>;
// these headers are mostly templates and inline class bodies, so they
// inherit the rename automatically.
#include "compositecodec.h"
#include "simdbitpacking.h"
#include "simdfastpfor.h"
#include "variablebyte.h"
#include "common.h"

#include <algorithm>
#include <cstring>
#include <vector>

#ifndef IKAFSSN_PFD_TIER_NAME
#error "IKAFSSN_PFD_TIER_NAME must be set per tier (sse42 / avx2 / avx512bw / avx512vbmi2)"
#endif

#define IKAFSSN_PFD_TIER_NS_(x) ikafssn_pfd_##x
#define IKAFSSN_PFD_TIER_NS(x)  IKAFSSN_PFD_TIER_NS_(x)

namespace ikafssn::pfd::IKAFSSN_PFD_TIER_NS(IKAFSSN_PFD_TIER_NAME) {

namespace {

namespace pfor_ns = FastPForLib;

constexpr int kBlockSize = 128;
constexpr int kBlockAlign = 16; // __m128i

using KixCodec = pfor_ns::CompositeCodec<
    pfor_ns::SIMDFastPFor<4>,
    pfor_ns::VariableByte>;

KixCodec& kix_codec() {
    thread_local KixCodec codec;
    return codec;
}

inline std::uint8_t bits_required(std::uint32_t v) {
    return v == 0 ? std::uint8_t(0) : std::uint8_t(32 - __builtin_clz(v));
}

// LSB-first variable bit-width packing for FOR-stream tails.  Each value
// is written into `bw` consecutive bits starting at the next bitpos in
// the output byte stream.  Total output: ceil(count*bw/8) bytes.
void pack_bits_lsb(const std::uint32_t* in, std::uint32_t count,
                   std::uint8_t bw, std::uint8_t* out) {
    if (bw == 0 || count == 0) return;
    const std::size_t total_bits  = std::size_t(count) * bw;
    const std::size_t total_bytes = (total_bits + 7) / 8;
    std::memset(out, 0, total_bytes);
    std::uint64_t bitpos = 0;
    for (std::uint32_t i = 0; i < count; i++) {
        std::uint64_t v = in[i];
        std::uint64_t bp = bitpos;
        std::uint8_t  bits_left = bw;
        while (bits_left > 0) {
            std::uint64_t byte_idx = bp >> 3;
            std::uint8_t  bit_idx  = std::uint8_t(bp & 7);
            std::uint8_t  bits_in_byte = std::uint8_t(8 - bit_idx);
            std::uint8_t  take = bits_left < bits_in_byte ? bits_left : bits_in_byte;
            std::uint64_t chunk_mask = (std::uint64_t(1) << take) - 1;
            std::uint8_t  chunk = std::uint8_t(v & chunk_mask);
            out[byte_idx] |= std::uint8_t(chunk << bit_idx);
            v >>= take;
            bp += take;
            bits_left = std::uint8_t(bits_left - take);
        }
        bitpos += bw;
    }
}

void unpack_bits_lsb(const std::uint8_t* in, std::uint32_t count,
                     std::uint8_t bw, std::uint32_t* out) {
    if (count == 0) return;
    if (bw == 0) {
        for (std::uint32_t i = 0; i < count; i++) out[i] = 0;
        return;
    }
    std::uint64_t bitpos = 0;
    for (std::uint32_t i = 0; i < count; i++) {
        std::uint64_t bp = bitpos;
        std::uint8_t  bits_left = bw;
        std::uint64_t v = 0;
        std::uint8_t  out_shift = 0;
        while (bits_left > 0) {
            std::uint64_t byte_idx = bp >> 3;
            std::uint8_t  bit_idx  = std::uint8_t(bp & 7);
            std::uint8_t  bits_in_byte = std::uint8_t(8 - bit_idx);
            std::uint8_t  take = bits_left < bits_in_byte ? bits_left : bits_in_byte;
            std::uint64_t chunk_mask = (std::uint64_t(1) << take) - 1;
            std::uint8_t  chunk = std::uint8_t((in[byte_idx] >> bit_idx) & std::uint8_t(chunk_mask));
            v |= std::uint64_t(chunk) << out_shift;
            out_shift = std::uint8_t(out_shift + take);
            bp += take;
            bits_left = std::uint8_t(bits_left - take);
        }
        out[i] = std::uint32_t(v);
        bitpos += bw;
    }
}

// FOR block: [u32 min][u8 b][3 B pad][16*b bytes bitpacked (value-min)].
// Total = 8 + 16*b bytes, body 8 B aligned within the block.
void encode_block_for(const std::uint32_t* in_128,
                         std::vector<std::uint8_t>& out) {
    std::uint32_t mn = in_128[0];
    std::uint32_t mx = in_128[0];
    for (int i = 1; i < kBlockSize; i++) {
        if (in_128[i] < mn) mn = in_128[i];
        if (in_128[i] > mx) mx = in_128[i];
    }
    const std::uint32_t spread = mx - mn;
    const std::uint8_t b = bits_required(spread);
    const std::size_t before = out.size();

    out.resize(before + 8);
    std::memcpy(out.data() + before, &mn, sizeof(std::uint32_t));
    out[before + 4] = b;
    out[before + 5] = 0;
    out[before + 6] = 0;
    out[before + 7] = 0;
    if (b == 0) return;

    alignas(kBlockAlign) std::uint32_t shifted[kBlockSize];
    for (int i = 0; i < kBlockSize; i++) shifted[i] = in_128[i] - mn;

    const std::size_t body_bytes = std::size_t(16) * b;
    out.resize(before + 8 + body_bytes);

    alignas(kBlockAlign) __m128i tmp[32];
    pfor_ns::simdpackwithoutmask(shifted, tmp, b);
    std::memcpy(out.data() + before + 8, tmp, body_bytes);
}

bool decode_block_for(const std::uint8_t*& p, const std::uint8_t* end,
                         std::uint32_t* out_128) {
    if (std::size_t(end - p) < 8) return false;
    std::uint32_t mn;
    std::memcpy(&mn, p, sizeof(std::uint32_t));
    const std::uint8_t b = p[4];
    if (b > 32) return false;
    const std::size_t body_bytes = std::size_t(16) * b;
    if (std::size_t(end - p) < 8 + body_bytes) return false;
    if (b == 0) {
        for (int i = 0; i < kBlockSize; i++) out_128[i] = mn;
        p += 8;
        return true;
    }
    alignas(kBlockAlign) __m128i tmp[32];
    std::memcpy(tmp, p + 8, body_bytes);
    pfor_ns::simdunpack(tmp, out_128, b);
    for (int i = 0; i < kBlockSize; i++) out_128[i] += mn;
    p += 8 + body_bytes;
    return true;
}

void encode_for_stream(const std::uint32_t* abs_pos_array,
                          std::uint32_t count,
                          std::vector<std::uint8_t>& out) {
    const std::uint32_t num_blocks = count / kBlockSize;
    const std::uint32_t tail_count = count % kBlockSize;

    for (std::uint32_t b = 0; b < num_blocks; b++) {
        encode_block_for(abs_pos_array + b * kBlockSize, out);
    }

    out.push_back(static_cast<std::uint8_t>(tail_count));
    if (tail_count == 0) return;

    const std::uint32_t* tail = abs_pos_array + num_blocks * kBlockSize;
    std::uint32_t mn = tail[0];
    std::uint32_t mx = tail[0];
    for (std::uint32_t i = 1; i < tail_count; i++) {
        if (tail[i] < mn) mn = tail[i];
        if (tail[i] > mx) mx = tail[i];
    }
    const std::uint8_t  tail_b     = bits_required(mx - mn);
    const std::size_t   body_bits  = std::size_t(tail_count) * tail_b;
    const std::size_t   body_bytes = (body_bits + 7) / 8;
    const std::size_t   before     = out.size();
    out.resize(before + 4 + 1 + body_bytes);
    std::memcpy(out.data() + before, &mn, sizeof(std::uint32_t));
    out[before + 4] = tail_b;
    if (body_bytes > 0) {
        // Pack (tail[i] - mn) values.
        std::uint32_t shifted[kBlockSize];
        for (std::uint32_t i = 0; i < tail_count; i++) shifted[i] = tail[i] - mn;
        pack_bits_lsb(shifted, tail_count, tail_b, out.data() + before + 5);
    }
}

bool decode_for_stream(const std::uint8_t*& p, const std::uint8_t* end,
                          std::uint32_t count, std::uint32_t* out) {
    const std::uint32_t num_blocks = count / kBlockSize;
    const std::uint32_t tail_count = count % kBlockSize;

    for (std::uint32_t b = 0; b < num_blocks; b++) {
        if (!decode_block_for(p, end, out + b * kBlockSize)) return false;
    }

    if (p >= end) return false;
    const std::uint8_t got_tail = *p++;
    if (got_tail != tail_count) return false;
    if (tail_count == 0) return true;

    if (std::size_t(end - p) < 5) return false;
    std::uint32_t tail_min;
    std::memcpy(&tail_min, p, sizeof(std::uint32_t));
    p += sizeof(std::uint32_t);
    const std::uint8_t tail_b = *p++;
    if (tail_b > 32) return false;
    const std::size_t body_bits  = std::size_t(tail_count) * tail_b;
    const std::size_t body_bytes = (body_bits + 7) / 8;
    if (std::size_t(end - p) < body_bytes) return false;

    std::uint32_t* tail_out = out + num_blocks * kBlockSize;
    unpack_bits_lsb(p, tail_count, tail_b, tail_out);
    for (std::uint32_t i = 0; i < tail_count; i++) tail_out[i] += tail_min;
    p += body_bytes;
    return true;
}

inline void set_kind_bits(std::uint8_t* kind_map, std::uint32_t i, std::uint8_t kind) {
    kind_map[i >> 2] |= std::uint8_t((kind & 0x03) << ((i & 3) * 2));
}

inline std::uint8_t get_kind_bits(const std::uint8_t* kind_map, std::uint32_t i) {
    return std::uint8_t((kind_map[i >> 2] >> ((i & 3) * 2)) & 0x03);
}

} // anonymous namespace (TU-private helpers)

// Count partition / short1 / short2 entries directly from the 2-bit
// kind map.  The kind_map encoding is:
//
//   00 = short_occ1     -> short1
//   01 = short_occ_ge2  -> short2
//   10 = partition      -> partition
//   11 = reserved
//
// Bulk path: process 32 pairs (8 bytes) per iteration via three
// popcount-of-mask operations, isolating each kind via bitwise
// expressions on the low and high bits of every pair.  Under -mavx2 /
// -mavx512bw the compiler auto-vectorises the AND/POPCNT chain into
// vpshufb-based byte popcount where applicable; on SSE4.2 the
// per-tier -mpopcnt produces hardware POPCNT for the inner u64.
void popcount_kinds(const std::uint8_t* km,
                    std::uint32_t distinct_count,
                    std::uint32_t* p_partition,
                    std::uint32_t* p_short1,
                    std::uint32_t* p_short2) noexcept {
    constexpr std::uint64_t kLow  = 0x5555555555555555ULL;
    std::uint32_t partition_count = 0;
    std::uint32_t short1_count    = 0;
    std::uint32_t short2_count    = 0;

    // Whole-chunk loop: 32 pairs (8 bytes) at a time, only when the
    // chunk lies entirely within distinct_count.
    std::uint32_t k = 0;
    while (k + 32 <= distinct_count) {
        std::uint64_t w;
        std::memcpy(&w, km + (k >> 2), 8);
        const std::uint64_t lo_bits = w & kLow;
        const std::uint64_t hi_bits = (w >> 1) & kLow;
        // 00 -> short1: ~hi & ~lo  (masked to pair positions via kLow)
        // 01 -> short2: ~hi &  lo
        // 10 -> partition: hi & ~lo
        // 11 -> reserved (caller must still invoke get_kind_bits to detect it)
        partition_count += static_cast<std::uint32_t>(
            __builtin_popcountll(hi_bits & ~lo_bits));
        short2_count    += static_cast<std::uint32_t>(
            __builtin_popcountll(~hi_bits & lo_bits & kLow));
        short1_count    += static_cast<std::uint32_t>(
            __builtin_popcountll(~hi_bits & ~lo_bits & kLow));
        k += 32;
    }

    // Tail: per-pair scalar.
    for (; k < distinct_count; ++k) {
        std::uint8_t kind = get_kind_bits(km, k);
        switch (kind) {
            case 0: short1_count++;    break;
            case 1: short2_count++;    break;
            case 2: partition_count++; break;
            default: break;  // 11 reserved; counted as "other" (caller validates)
        }
    }

    *p_partition = partition_count;
    *p_short1    = short1_count;
    *p_short2    = short2_count;
}

// ===== .kix encode (distinct seq_id delta stream → SIMDFastPFor + VByte tail) =====
//
// Only the leading `[u32 distinct_count]` is written; the body length
// is derived at decode time from the EF dictionary's posting byte
// length.  On-disk per-posting-list layout:
//   [u32 distinct_count]   — number of distinct seq_ids represented
//   [u32 body[N]]          — codec output (N = (bytes - 4) / 4)

std::size_t encode_posting_kix(const std::uint32_t* delta_array,
                               std::uint32_t count,
                               std::vector<std::uint8_t>& out) {
    const std::size_t before = out.size();
    out.resize(before + 4);
    std::memcpy(out.data() + before, &count, sizeof(std::uint32_t));
    if (count == 0) {
        return 4;
    }

    const std::size_t worst_words = std::size_t(count) * 2 + 1024;
    std::vector<std::uint32_t> codec_out(worst_words);
    std::size_t nvalue = codec_out.size();
    kix_codec().encodeArray(delta_array, count, codec_out.data(), nvalue);

    const std::size_t body_bytes = nvalue * sizeof(std::uint32_t);
    out.resize(out.size() + body_bytes);
    std::memcpy(out.data() + before + 4, codec_out.data(), body_bytes);

    return out.size() - before;
}

// ===== .kpx encode =====
//
// The body starts directly at the 2-bit kind map; counts are derived
// at decode time from the kind map (popcount_kinds for partition /
// short1 / short2) and by summing the u8 occ_count[] array
// (short2_position_count).  See src/index/pfd_codec.hpp for the wire
// format.
//
// Empty posting lists (distinct_count == 0) emit zero bytes; the
// caller's dictionary is responsible for delimiting each posting list
// (an empty .kpx posting list aliases the next non-empty posting
// list's start).

std::size_t encode_posting_kpx(const std::uint32_t* /*distinct_sid*/,
                               const std::uint32_t* occ_count,
                               std::uint32_t distinct_count,
                               const std::uint32_t* abs_pos_array,
                               std::uint32_t /*position_count*/,
                               std::uint32_t freq_threshold_part,
                               std::vector<std::uint8_t>& out) {
    const std::size_t before = out.size();

    if (distinct_count == 0) {
        // Empty posting list emits 0 bytes (no header).  The decoder
        // gates on kix_count > 0 and never reads this region.
        return 0;
    }

    // Pass 1: classify each distinct sid by occ_count.
    std::vector<std::uint8_t> kinds(distinct_count);
    std::uint32_t partition_count = 0, short1_count = 0, short2_count = 0;
    std::uint32_t short2_position_count = 0;
    for (std::uint32_t k = 0; k < distinct_count; k++) {
        const std::uint32_t occ = occ_count[k];
        std::uint8_t kind;
        if (occ > freq_threshold_part) {
            kind = 2;  // partition
            partition_count++;
        } else if (occ == 1) {
            kind = 0;  // short_occ1
            short1_count++;
        } else {
            kind = 1;  // short_occ_ge2
            short2_count++;
            short2_position_count += occ;
        }
        kinds[k] = kind;
    }

    // 2-bit kind map (the body now starts here — no preceding u32 fields).
    const std::size_t kind_map_bytes = (std::size_t(distinct_count) * 2 + 7) / 8;
    const std::size_t kind_map_off = out.size();
    out.resize(out.size() + kind_map_bytes);
    std::memset(out.data() + kind_map_off, 0, kind_map_bytes);
    for (std::uint32_t k = 0; k < distinct_count; k++) {
        set_kind_bits(out.data() + kind_map_off, k, kinds[k]);
    }

    // Pass 2: emit partition groups in distinct_sid order.
    std::uint32_t pos_cursor = 0;
    for (std::uint32_t k = 0; k < distinct_count; k++) {
        const std::uint32_t cnt = occ_count[k];
        if (kinds[k] == 2) {
            const std::size_t off = out.size();
            out.resize(off + sizeof(std::uint32_t));
            std::memcpy(out.data() + off, &cnt, sizeof(std::uint32_t));
            encode_for_stream(abs_pos_array + pos_cursor, cnt, out);
        }
        pos_cursor += cnt;
    }

    // Pass 3: short_occ1 — concatenated 1-position-per-cluster FOR stream.
    if (short1_count > 0) {
        std::vector<std::uint32_t> short1_buf;
        short1_buf.reserve(short1_count);
        std::uint32_t cursor = 0;
        for (std::uint32_t k = 0; k < distinct_count; k++) {
            if (kinds[k] == 0) {
                short1_buf.push_back(abs_pos_array[cursor]);
            }
            cursor += occ_count[k];
        }
        encode_for_stream(short1_buf.data(), short1_count, out);
    }

    // Pass 4: short_occ_ge2 — u8 occ_count[] + concatenated FOR stream.
    // occ_count[k] is in [2, freq_threshold_part] (<= 255) for short_occ_ge2
    // entries, so the static_cast<uint8_t> never loses information.
    if (short2_count > 0) {
        const std::size_t occ_off = out.size();
        out.resize(occ_off + short2_count);
        std::uint32_t out_idx = 0;
        for (std::uint32_t k = 0; k < distinct_count; k++) {
            if (kinds[k] == 1) {
                out[occ_off + out_idx++] = static_cast<std::uint8_t>(occ_count[k]);
            }
        }

        std::vector<std::uint32_t> short2_buf;
        short2_buf.reserve(short2_position_count);
        std::uint32_t cursor = 0;
        for (std::uint32_t k = 0; k < distinct_count; k++) {
            if (kinds[k] == 1) {
                for (std::uint32_t e = 0; e < occ_count[k]; e++) {
                    short2_buf.push_back(abs_pos_array[cursor + e]);
                }
            }
            cursor += occ_count[k];
        }
        encode_for_stream(short2_buf.data(), short2_position_count, out);
    }

    return out.size() - before;
}

// ===== open_stream_kix: decode the entire .kix posting list into the StreamCtx =====

bool open_stream_kix(const std::uint8_t* posting_list, std::size_t bytes,
                     ikafssn::pfd::StreamCtx& ctx) {
    ctx.count = 0;
    if (bytes == 0) return true;
    if (bytes < 4) return false;  // header is just [u32 distinct_count]

    std::uint32_t count;
    std::memcpy(&count, posting_list, sizeof(std::uint32_t));
    if (count == 0) return true;

    // body_words is derived from the posting list byte length supplied
    // by the EF dictionary (no on-wire body_words).
    const std::size_t body_bytes_avail = bytes - 4;
    if ((body_bytes_avail % sizeof(std::uint32_t)) != 0) return false;
    const std::uint32_t body_words = static_cast<std::uint32_t>(
        body_bytes_avail / sizeof(std::uint32_t));

    // Encoder writes posting lists as [u32 header][u32-aligned body], and
    // each posting list is itself an integer number of u32 words, so the
    // body offset relative to the posting file base is always 4-byte
    // aligned. When the underlying pointer is u32-aligned we hand it to
    // FastPFor directly to skip the per-call memcpy + heap allocation.
    // Otherwise we fall back to a per-thread scratch vector that grows
    // monotonically — never per-call heap allocation.
    const std::uint8_t* body_bytes_ptr = posting_list + 4;
    const std::uint32_t* codec_in_ptr;
    if ((reinterpret_cast<std::uintptr_t>(body_bytes_ptr) % alignof(std::uint32_t)) == 0) {
        codec_in_ptr = reinterpret_cast<const std::uint32_t*>(body_bytes_ptr);
    } else {
        thread_local std::vector<std::uint32_t> codec_in_scratch;
        if (codec_in_scratch.size() < body_words) codec_in_scratch.resize(body_words);
        std::memcpy(codec_in_scratch.data(), body_bytes_ptr, body_bytes_avail);
        codec_in_ptr = codec_in_scratch.data();
    }

    // `decoded` grows to a high-water mark: its size() is never read — the
    // valid range is [0, ctx.count) — so a shorter posting list reuses the
    // buffer as it stands and the decode below overwrites what it needs.
    if (ctx.decoded.size() < count) ctx.decoded.resize(count);
    std::size_t nvalue = count;
    kix_codec().decodeArray(codec_in_ptr, body_words,
                            ctx.decoded.data(), nvalue);
    if (nvalue != count) return false;

    // Encoder writes [abs_first, d1, d2, ...] over the **distinct**
    // seq_id stream; cumulative sum reconstructs absolute distinct seq_ids.
    for (std::uint32_t i = 1; i < count; i++) {
        ctx.decoded[i] += ctx.decoded[i - 1];
    }

    ctx.count = count;
    return true;
}

// ===== open_stream_kpx_for_candidates: candidate-set-driven decode =====

namespace {

using Run      = ikafssn::pfd::PosDecodeScratch::Run;
using Selected = ikafssn::pfd::PosDecodeScratch::Selected;

// The range popcount's bulk loop needs 32 byte-aligned entries before it
// does any work, so a shorter gap only pays its prologue and the call.
constexpr std::uint32_t kRangePopcountMinGap = 32;

// First index in [from, n) whose value is >= target, located by exponential
// probing followed by a binary search over the bracketed range.  Returns n
// when every remaining value is below target.
inline std::size_t gallop_lower_bound(const std::uint32_t* a, std::size_t from,
                                      std::size_t n, std::uint32_t target) {
    if (from >= n) return n;
    if (a[from] >= target) return from;
    std::size_t lo = from;          // invariant: a[lo] < target
    std::size_t step = 1;
    while (lo + step < n && a[lo + step] < target) {
        lo += step;
        step <<= 1;
    }
    const std::size_t hi = std::min(lo + step + 1, n);
    return std::size_t(std::lower_bound(a + lo + 1, a + hi, target) - a);
}

// Bring the per-kind ranks up to date over the kind map entries [first, last)
// without visiting them one at a time.
inline void advance_ranks(const std::uint8_t* km,
                          std::uint32_t first, std::uint32_t last,
                          std::uint32_t& r_partition,
                          std::uint32_t& r_short1,
                          std::uint32_t& r_short2) {
    std::uint32_t k = first;
    if (last - first >= kRangePopcountMinGap) {
        for (; k < last && (k & 3) != 0; k++) {
            switch (get_kind_bits(km, k)) {
                case 0: r_short1++;    break;
                case 1: r_short2++;    break;
                case 2: r_partition++; break;
                default: break;
            }
        }
        std::uint32_t np, n1, n2;
        popcount_kinds(km + (k >> 2), last - k, &np, &n1, &n2);
        r_partition += np;
        r_short1    += n1;
        r_short2    += n2;
        return;
    }
    for (; k < last; k++) {
        switch (get_kind_bits(km, k)) {
            case 0: r_short1++;    break;
            case 1: r_short2++;    break;
            case 2: r_partition++; break;
            default: break;
        }
    }
}

// Walk past one FOR block through its self-describing header.
inline bool skip_block_for(const std::uint8_t*& p, const std::uint8_t* end) {
    if (std::size_t(end - p) < 8) return false;
    const std::uint8_t b = p[4];
    if (b > 32) return false;
    const std::size_t body_bytes = std::size_t(16) * b;
    if (std::size_t(end - p) < 8 + body_bytes) return false;
    p += 8 + body_bytes;
    return true;
}

// Walk past a whole FOR stream of `count` elements without decoding any of it.
bool skip_for_stream(const std::uint8_t*& p, const std::uint8_t* end,
                     std::uint32_t count, std::uint64_t& skipped_blocks) {
    const std::uint32_t num_blocks = count / kBlockSize;
    const std::uint32_t tail_count = count % kBlockSize;

    for (std::uint32_t b = 0; b < num_blocks; b++) {
        if (!skip_block_for(p, end)) return false;
        skipped_blocks++;
    }

    if (p >= end) return false;
    if (*p++ != tail_count) return false;
    if (tail_count == 0) return true;

    if (std::size_t(end - p) < 5) return false;
    const std::uint8_t tail_b = p[4];
    if (tail_b > 32) return false;
    const std::size_t body_bytes = (std::size_t(tail_count) * tail_b + 7) / 8;
    if (std::size_t(end - p) < 5 + body_bytes) return false;
    p += 5 + body_bytes;
    skipped_blocks++;
    return true;
}

// Decode only the elements covered by `runs` out of a FOR stream of `count`
// elements, writing them to `out` in stream order.  `runs` must be ascending,
// disjoint and inside [0, count), and `out` must hold their total length.
// Blocks holding no selected element are walked past through their header.
bool decode_for_stream_selected(const std::uint8_t*& p, const std::uint8_t* end,
                                std::uint32_t count,
                                const Run* runs, std::size_t n_runs,
                                std::uint32_t* out,
                                std::uint64_t& skipped_blocks) {
    const std::uint32_t num_blocks = count / kBlockSize;
    const std::uint32_t tail_count = count % kBlockSize;

    std::size_t ri = 0;             // run under consideration
    std::uint32_t consumed = 0;     // elements of runs[ri] already emitted
    alignas(kBlockAlign) std::uint32_t buf[kBlockSize];

    // Drop the runs that end before this block starts.
    auto seek_runs = [&](std::uint32_t base) {
        while (ri < n_runs && runs[ri].first + runs[ri].length <= base) {
            ri++;
            consumed = 0;
        }
    };
    // Emit every selected element of the decoded block spanning [base, lim).
    auto emit = [&](std::uint32_t base, std::uint32_t lim) {
        while (ri < n_runs && runs[ri].first < lim) {
            const std::uint32_t run_end = runs[ri].first + runs[ri].length;
            const std::uint32_t hi = run_end < lim ? run_end : lim;
            for (std::uint32_t e = runs[ri].first + consumed; e < hi; e++) {
                *out++ = buf[e - base];
            }
            if (run_end > lim) {          // continues into the next block
                consumed = lim - runs[ri].first;
                return;
            }
            ri++;
            consumed = 0;
        }
    };

    for (std::uint32_t b = 0; b < num_blocks; b++) {
        const std::uint32_t base = b * kBlockSize;
        const std::uint32_t lim  = base + kBlockSize;
        seek_runs(base);
        if (ri >= n_runs || runs[ri].first >= lim) {
            if (!skip_block_for(p, end)) return false;
            skipped_blocks++;
            continue;
        }
        if (!decode_block_for(p, end, buf)) return false;
        emit(base, lim);
    }

    if (p >= end) return false;
    if (*p++ != tail_count) return false;
    if (tail_count == 0) return true;

    if (std::size_t(end - p) < 5) return false;
    std::uint32_t tail_min;
    std::memcpy(&tail_min, p, sizeof(std::uint32_t));
    const std::uint8_t tail_b = p[4];
    if (tail_b > 32) return false;
    const std::size_t body_bytes = (std::size_t(tail_count) * tail_b + 7) / 8;
    if (std::size_t(end - p) < 5 + body_bytes) return false;

    const std::uint32_t base = num_blocks * kBlockSize;
    const std::uint32_t lim  = base + tail_count;
    seek_runs(base);
    if (ri >= n_runs || runs[ri].first >= lim) {
        p += 5 + body_bytes;
        skipped_blocks++;
        return true;
    }
    unpack_bits_lsb(p + 5, tail_count, tail_b, buf);
    for (std::uint32_t i = 0; i < tail_count; i++) buf[i] += tail_min;
    p += 5 + body_bytes;
    emit(base, lim);
    return true;
}

} // anonymous namespace (selection-pass helpers)

bool open_stream_kpx_for_candidates(
        const std::uint8_t* posting_list, std::size_t bytes,
        const std::uint32_t* kix_decoded, std::size_t kix_count,
        const std::uint32_t* candidates, std::size_t n_candidates,
        ikafssn::pfd::PosDecodeScratch& scratch) {

    auto& out_cand = scratch.out_candidate_idx;
    auto& out_off  = scratch.out_offsets;
    auto& out_pos  = scratch.out_positions;
    out_cand.clear();
    out_pos.clear();
    out_off.assign(1, 0);
    scratch.skipped_kind_entries = 0;
    scratch.skipped_partition_groups = 0;
    scratch.skipped_for_blocks = 0;

    if (n_candidates == 0) return true;
    if (bytes == 0) return true;
    const std::uint32_t distinct_count = static_cast<std::uint32_t>(kix_count);
    if (distinct_count == 0) return true;

    const std::uint8_t* p = posting_list;
    const std::uint8_t* const end = posting_list + bytes;

    const std::size_t kind_map_bytes = (std::size_t(distinct_count) * 2 + 7) / 8;
    if (std::size_t(end - p) < kind_map_bytes) return false;
    const std::uint8_t* const kind_map = p;
    p += kind_map_bytes;

    std::uint32_t partition_count, short1_count, short2_count;
    popcount_kinds(kind_map, distinct_count,
                   &partition_count, &short1_count, &short2_count);
    // A reserved (11) entry is counted into none of the three, so this one
    // comparison stands in for a per-entry validity check.  It is what lets
    // the selection pass below treat "neither partition nor short_occ1" as
    // short_occ_ge2 without re-testing every entry.
    if (distinct_count != partition_count + short1_count + short2_count) return false;

    // Selection pass: galloping intersection of the .kix distinct sid array
    // with the candidate array.  The lagging side jumps ahead by exponential
    // probing, so the kind map is resolved one entry at a time only where
    // the two meet.
    auto& selected = scratch.selected;
    selected.clear();
    std::size_t i = 0, ci = 0;
    std::uint32_t prev_i = 0;   // ranks are up to date over [0, prev_i)
    std::uint32_t r_partition = 0, r_short1 = 0, r_short2 = 0;
    std::uint32_t n_sel_short1 = 0, n_sel_short2 = 0;
    while (i < distinct_count && ci < n_candidates) {
        const std::uint32_t sid  = kix_decoded[i];
        const std::uint32_t cand = candidates[ci];
        if (sid < cand) {
            i = gallop_lower_bound(kix_decoded, i + 1, distinct_count, cand);
            continue;
        }
        if (cand < sid) {
            ci = gallop_lower_bound(candidates, ci + 1, n_candidates, sid);
            continue;
        }
        const std::uint32_t idx = static_cast<std::uint32_t>(i);
        scratch.skipped_kind_entries += idx - prev_i;
        advance_ranks(kind_map, prev_i, idx, r_partition, r_short1, r_short2);
        const std::uint8_t kind = get_kind_bits(kind_map, idx);
        std::uint32_t rank;
        if (kind == 2)      { rank = r_partition++; }
        else if (kind == 0) { rank = r_short1++; n_sel_short1++; }
        else                { rank = r_short2++; n_sel_short2++; }
        selected.push_back(Selected{static_cast<std::uint32_t>(ci), rank, kind});
        prev_i = idx + 1;
        i++;
        ci++;
    }
    scratch.skipped_kind_entries += distinct_count - prev_i;

    // Nothing here can produce output, so no stream is touched at all.
    if (selected.empty()) return true;

    // Partition groups.  Every group header has to be read to find the next
    // one, but a group no candidate asked for has its FOR stream walked past.
    auto& part_pos = scratch.selected_partition_positions;
    auto& part_off = scratch.selected_partition_offsets;
    part_pos.clear();
    part_off.assign(1, 0);
    if (partition_count > 0) {
        std::size_t si = 0;
        while (si < selected.size() && selected[si].kind != 2) si++;
        std::uint32_t next_rank =
            (si < selected.size()) ? selected[si].rank : UINT32_MAX;
        for (std::uint32_t g = 0; g < partition_count; g++) {
            if (std::size_t(end - p) < sizeof(std::uint32_t)) return false;
            std::uint32_t gcnt;
            std::memcpy(&gcnt, p, sizeof(std::uint32_t));
            p += sizeof(std::uint32_t);
            if (gcnt == 0) return false;
            if (g != next_rank) {
                if (!skip_for_stream(p, end, gcnt, scratch.skipped_for_blocks)) return false;
                scratch.skipped_partition_groups++;
                continue;
            }
            const std::size_t base = part_pos.size();
            part_pos.resize(base + gcnt);
            if (!decode_for_stream(p, end, gcnt, part_pos.data() + base)) return false;
            part_off.push_back(static_cast<std::uint32_t>(part_pos.size()));
            si++;
            while (si < selected.size() && selected[si].kind != 2) si++;
            next_rank = (si < selected.size()) ? selected[si].rank : UINT32_MAX;
        }
    }

    // short_occ1 sub-bucket: one position per entry, so a selected rank is
    // its own single-element run.
    auto& short1_pos = scratch.selected_short1_positions;
    auto& runs = scratch.selected_runs;
    if (short1_count > 0 && n_sel_short1 > 0) {
        // High-water mark, so a shorter posting list never re-zeroes it.
        if (short1_pos.size() < n_sel_short1) short1_pos.resize(n_sel_short1);
        if (n_sel_short1 == short1_count) {
            // Everything is wanted: running whole blocks straight through
            // beats selecting inside them, and needs no run list.
            if (!decode_for_stream(p, end, short1_count, short1_pos.data())) return false;
        } else {
            runs.clear();
            for (const auto& s : selected) {
                if (s.kind == 0) runs.push_back(Run{s.rank, 1});
            }
            if (!decode_for_stream_selected(p, end, short1_count,
                                            runs.data(), runs.size(),
                                            short1_pos.data(),
                                            scratch.skipped_for_blocks)) {
                return false;
            }
        }
    } else if (short1_count > 0 && n_sel_short2 > 0) {
        // Nothing wanted here, but the stream sits between the partition
        // groups and the short_occ_ge2 sub-bucket.
        if (!skip_for_stream(p, end, short1_count, scratch.skipped_for_blocks)) return false;
    }

    // short_occ_ge2 sub-bucket: u8 occ_count[] then the FOR stream.  The
    // occ_count[] array is read in full because a selected entry's run starts
    // at the prefix sum of every entry before it.  Being the last region, it
    // is skipped outright when no candidate selected one of its entries.
    auto& short2_pos = scratch.selected_short2_positions;
    auto& short2_occ = scratch.short2_occ;
    if (n_sel_short2 > 0) {
        if (std::size_t(end - p) < short2_count) return false;
        short2_occ.resize(short2_count);
        std::memcpy(short2_occ.data(), p, short2_count);
        p += short2_count;

        runs.clear();
        std::size_t si = 0;
        while (si < selected.size() && selected[si].kind != 1) si++;
        std::uint32_t cum = 0;
        std::uint32_t sel_positions = 0;
        for (std::uint32_t r = 0; r < short2_count; r++) {
            const std::uint8_t occ = short2_occ[r];
            if (occ < 2) return false;
            if (si < selected.size() && selected[si].rank == r) {
                runs.push_back(Run{cum, occ});
                sel_positions += occ;
                si++;
                while (si < selected.size() && selected[si].kind != 1) si++;
            }
            cum += occ;
        }
        if (short2_pos.size() < sel_positions) short2_pos.resize(sel_positions);
        if (!decode_for_stream_selected(p, end, cum, runs.data(), runs.size(),
                                        short2_pos.data(), scratch.skipped_for_blocks)) {
            return false;
        }
    }

    // Assemble the sparse CSR.  Each kind's selected positions landed in
    // selection order, so one cursor per kind walks them in step with the
    // ascending candidate-index order of `selected`.
    out_cand.reserve(selected.size());
    out_off.reserve(selected.size() + 1);
    std::size_t c_partition = 0, c_short1 = 0, c_short2 = 0;
    for (const auto& s : selected) {
        out_cand.push_back(s.candidate_idx);
        if (s.kind == 2) {
            const std::uint32_t lo = part_off[c_partition];
            const std::uint32_t hi = part_off[c_partition + 1];
            out_pos.insert(out_pos.end(),
                           part_pos.begin() + lo, part_pos.begin() + hi);
            c_partition++;
        } else if (s.kind == 0) {
            out_pos.push_back(short1_pos[c_short1++]);
        } else {
            const std::uint32_t occ = short2_occ[s.rank];
            out_pos.insert(out_pos.end(),
                           short2_pos.begin() + c_short2,
                           short2_pos.begin() + c_short2 + occ);
            c_short2 += occ;
        }
        out_off.push_back(static_cast<std::uint32_t>(out_pos.size()));
    }
    return true;
}

} // namespace ikafssn::pfd::ikafssn_pfd_<tier>
