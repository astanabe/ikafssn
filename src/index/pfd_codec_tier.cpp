// Phase 6 — per-tier custom block codec.
//
// This translation unit is compiled FOUR TIMES on x86_64 (once per ISA
// tier) and once on aarch64 via the ikafssn_pfd_<tier> OBJECT libraries
// declared in the top-level CMakeLists.txt.  Each compilation is given:
//
//   -DFastPForLib=FastPForLib_<tier>     (renames FastPFor's namespace
//                                          at preprocessor time so the
//                                          per-tier sets of symbols do
//                                          not collide at link time)
//   -DIKAFSSN_PFD_TIER_NAME=<tier>       (used here to name the
//                                          tier-specific ikafssn::pfd
//                                          inner namespace)
//   -m<arch> ...                          (controls the instructions the
//                                          bitpacker primitives and the
//                                          surrounding scalar code in
//                                          this file are allowed to use)
//
// For .kix posting lists (sorted distinct seq_id delta stream — wire format
// unchanged from v7) we drive FastPFor's CompositeCodec<SIMDFastPFor<4>,
// VariableByte> directly.
//
// For .kpx posting lists (absolute position stream) Phase 6 keeps the per-
// (kmer, seq_id) partition grouping but classifies every distinct
// seq_id via a 2-bit kind map instead of carrying redundant seq_id
// bytes.  The short bucket is also split into occ=1 and occ>=2 sub-
// buckets so single-position clusters drop the u8 occ_count entirely.
// FOR-block headers are widened to 8 B (proposal F) and stream tails
// switch from a varint stream to a packed bit-width stream (proposal D).

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

// v8 FOR block: [u32 min][u8 b][3 B pad][16*b bytes bitpacked (value-min)].
// Total = 8 + 16*b bytes, body 8 B aligned within the block.
void encode_block_for_v8(const std::uint32_t* in_128,
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

bool decode_block_for_v8(const std::uint8_t*& p, const std::uint8_t* end,
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

void encode_for_stream_v8(const std::uint32_t* abs_pos_array,
                          std::uint32_t count,
                          std::vector<std::uint8_t>& out) {
    const std::uint32_t num_blocks = count / kBlockSize;
    const std::uint32_t tail_count = count % kBlockSize;

    for (std::uint32_t b = 0; b < num_blocks; b++) {
        encode_block_for_v8(abs_pos_array + b * kBlockSize, out);
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

bool decode_for_stream_v8(const std::uint8_t*& p, const std::uint8_t* end,
                          std::uint32_t count, std::uint32_t* out) {
    const std::uint32_t num_blocks = count / kBlockSize;
    const std::uint32_t tail_count = count % kBlockSize;

    for (std::uint32_t b = 0; b < num_blocks; b++) {
        if (!decode_block_for_v8(p, end, out + b * kBlockSize)) return false;
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

// Phase 7d dedup C: count partition / short1 / short2 entries directly
// from the 2-bit kind map.  Replaces the redundant
// [u32 partition_count][u32 short1_count][u32 short2_count] header
// fields.  The kind_map encoding is:
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

// Phase 7d dedup D: sum a u8 array — replaces the redundant
// [u32 short2_position_count] header field.  Compiler auto-vectorises
// to vpsadbw under AVX2 / AVX512BW.  Tier-namespaced (not anonymous)
// for the same VTable-dispatch reason as popcount_kinds.
std::uint32_t horizontal_sum_u8(const std::uint8_t* arr,
                                std::uint32_t n) noexcept {
    std::uint32_t sum = 0;
    for (std::uint32_t i = 0; i < n; ++i) sum += arr[i];
    return sum;
}

// ===== .kix encode (distinct seq_id delta stream → SIMDFastPFor + VByte tail) =====
//
// Phase 7c dedup B: only the leading `[u32 distinct_count]` is written;
// `body_words` is derived at decode time from the EF dictionary's
// posting_byte_length.  On-disk per-posting-list layout:
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
// Phase 7c+7d: all four redundant header fields removed.  The body
// starts directly at the 2-bit kind map; counts are derived at decode
// time from the kind map (popcount_kinds for partition/short1/short2)
// and from the u8 occ_count[] array (horizontal_sum_u8 for
// short2_position_count).  See src/index/pfd_codec.hpp for the wire
// format.
//
// Empty posting lists (distinct_count == 0) emit zero bytes; the
// caller's offset table is responsible for delimiting the per-k-mer
// region (an empty .kpx posting list aliases the next non-empty
// posting list's start).

std::size_t encode_posting_kpx(const std::uint32_t* /*distinct_sid*/,
                               const std::uint32_t* occ_count,
                               std::uint32_t distinct_count,
                               const std::uint32_t* abs_pos_array,
                               std::uint32_t /*position_count*/,
                               std::uint32_t freq_threshold_part,
                               std::vector<std::uint8_t>& out) {
    const std::size_t before = out.size();

    if (distinct_count == 0) {
        // Phase 7d: empty posting list emits 0 bytes (no header).  The
        // decoder gates on kix_count > 0 and never reads this region.
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
            encode_for_stream_v8(abs_pos_array + pos_cursor, cnt, out);
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
        encode_for_stream_v8(short1_buf.data(), short1_count, out);
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
        encode_for_stream_v8(short2_buf.data(), short2_position_count, out);
    }

    return out.size() - before;
}

// ===== open_stream_kix: decode the entire .kix posting list list into the StreamCtx =====

bool open_stream_kix(const std::uint8_t* posting_list, std::size_t bytes,
                     ikafssn::pfd::StreamCtx& ctx) {
    ctx.decoded.clear();
    ctx.count = 0;
    ctx.pos = 0;
    if (bytes == 0) return true;
    if (bytes < 4) return false;  // Phase 7c: header is just [u32 distinct_count]

    std::uint32_t count;
    std::memcpy(&count, posting_list, sizeof(std::uint32_t));
    if (count == 0) return true;

    // Phase 7c dedup B: body_words derived from the posting list byte
    // length supplied by the EF dictionary (no on-wire body_words).
    const std::size_t body_bytes_avail = bytes - 4;
    if ((body_bytes_avail % sizeof(std::uint32_t)) != 0) return false;
    const std::uint32_t body_words = static_cast<std::uint32_t>(
        body_bytes_avail / sizeof(std::uint32_t));

    std::vector<std::uint32_t> codec_in(body_words);
    std::memcpy(codec_in.data(), posting_list + 4, body_bytes_avail);

    ctx.decoded.resize(count);
    std::size_t nvalue = ctx.decoded.size();
    kix_codec().decodeArray(codec_in.data(), body_words,
                            ctx.decoded.data(), nvalue);
    if (nvalue != count) return false;

    // Encoder writes [abs_first, d1, d2, ...] over the **distinct**
    // seq_id stream; cumulative sum reconstructs absolute distinct seq_ids.
    for (std::uint32_t i = 1; i < count; i++) {
        ctx.decoded[i] += ctx.decoded[i - 1];
    }

    ctx.count = count;
    ctx.pos = 0;
    return true;
}

// ===== open_stream_kpx_for_candidates: candidate-set-driven decode (v8) =====

bool open_stream_kpx_for_candidates(
        const std::uint8_t* posting_list, std::size_t bytes,
        const std::uint32_t* kix_decoded, std::size_t kix_count,
        const std::uint32_t* candidates, std::size_t n_candidates,
        ikafssn::pfd::PosDecodeScratch& scratch,
        std::vector<std::vector<std::uint32_t>>& out) {

    out.assign(n_candidates, std::vector<std::uint32_t>{});

    if (n_candidates == 0) return true;
    if (bytes == 0) return true;
    // Phase 7c+7d: all four redundant header fields removed.  Body
    // starts directly at the 2-bit kind map; per-kind counts are
    // derived via popcount_kinds, and short2_position_count via
    // horizontal_sum_u8 over the u8 occ_count[] array.
    const std::uint32_t distinct_count = static_cast<std::uint32_t>(kix_count);
    if (distinct_count == 0) return true;

    const std::uint8_t* p = posting_list;
    const std::uint8_t* end = posting_list + bytes;

    const std::size_t kind_map_bytes = (std::size_t(distinct_count) * 2 + 7) / 8;
    if (std::size_t(end - p) < kind_map_bytes) return false;
    const std::uint8_t* kind_map = p;
    p += kind_map_bytes;

    std::uint32_t partition_count, short1_count, short2_count;
    popcount_kinds(kind_map, distinct_count,
                   &partition_count, &short1_count, &short2_count);
    if (distinct_count != partition_count + short1_count + short2_count) return false;

    // Decode partition groups into scratch.partition_positions, building
    // partition_offsets[] alongside.
    auto& part_pos = scratch.partition_positions;
    auto& part_off = scratch.partition_offsets;
    part_pos.clear();
    part_off.assign(std::size_t(partition_count) + 1, 0);
    for (std::uint32_t g = 0; g < partition_count; g++) {
        if (std::size_t(end - p) < sizeof(std::uint32_t)) return false;
        std::uint32_t gcnt;
        std::memcpy(&gcnt, p, sizeof(std::uint32_t));
        p += sizeof(std::uint32_t);
        if (gcnt == 0) return false;
        part_off[g] = static_cast<std::uint32_t>(part_pos.size());
        const std::size_t base = part_pos.size();
        part_pos.resize(base + gcnt);
        if (!decode_for_stream_v8(p, end, gcnt, part_pos.data() + base)) return false;
    }
    part_off[partition_count] = static_cast<std::uint32_t>(part_pos.size());

    // Decode short_occ1 sub-bucket.
    auto& short1_pos = scratch.short1_positions;
    short1_pos.assign(short1_count, 0);
    if (short1_count > 0) {
        if (!decode_for_stream_v8(p, end, short1_count, short1_pos.data())) return false;
    }

    // Decode short_occ_ge2 sub-bucket: u8 occ_count[] then FOR stream.
    // Phase 7d dedup D: short2_position_count is derived as a horizontal
    // sum of the u8 occ_count[] array (the v8 in-blob u32 was redundant).
    auto& short2_occ = scratch.short2_occ;
    auto& short2_off = scratch.short2_offsets;
    auto& short2_pos = scratch.short2_positions;
    short2_occ.assign(short2_count, 0);
    short2_off.assign(std::size_t(short2_count) + 1, 0);
    std::uint32_t short2_position_count = 0;
    if (short2_count > 0) {
        if (std::size_t(end - p) < short2_count) return false;
        std::memcpy(short2_occ.data(), p, short2_count);
        p += short2_count;
        std::uint32_t cum = 0;
        for (std::uint32_t i = 0; i < short2_count; i++) {
            if (short2_occ[i] < 2) return false;
            short2_off[i] = cum;
            cum += short2_occ[i];
        }
        short2_off[short2_count] = cum;
        short2_position_count = cum;
    }
    short2_pos.assign(short2_position_count, 0);
    if (short2_position_count > 0) {
        if (!decode_for_stream_v8(p, end, short2_position_count, short2_pos.data())) {
            return false;
        }
    }

    // 2-pointer merge walk over kix_decoded × candidates, ranking each
    // distinct sid by its kind to pull positions out of the per-kind
    // decoded buffers.
    std::size_t ci = 0;
    std::uint32_t r_part = 0, r_short1 = 0, r_short2 = 0;
    for (std::uint32_t i = 0; i < distinct_count; i++) {
        const std::uint8_t kind = get_kind_bits(kind_map, i);
        if (kind == 3) return false;  // reserved
        const std::uint32_t sid = kix_decoded[i];
        while (ci < n_candidates && candidates[ci] < sid) ci++;
        const bool match = (ci < n_candidates && candidates[ci] == sid);
        if (kind == 2) {
            if (match) {
                const std::uint32_t lo = part_off[r_part];
                const std::uint32_t hi = part_off[r_part + 1];
                out[ci].assign(part_pos.begin() + lo, part_pos.begin() + hi);
            }
            r_part++;
        } else if (kind == 0) {
            if (match) {
                out[ci].assign(1, short1_pos[r_short1]);
            }
            r_short1++;
        } else { // kind == 1
            if (match) {
                const std::uint32_t lo = short2_off[r_short2];
                const std::uint32_t hi = short2_off[r_short2 + 1];
                out[ci].assign(short2_pos.begin() + lo, short2_pos.begin() + hi);
            }
            r_short2++;
        }
    }
    if (r_part   != partition_count) return false;
    if (r_short1 != short1_count)    return false;
    if (r_short2 != short2_count)    return false;
    return true;
}

} // namespace ikafssn::pfd::ikafssn_pfd_<tier>
