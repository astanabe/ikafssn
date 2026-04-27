// Phase 5e — per-tier custom block codec.
//
// This translation unit is compiled THREE TIMES (once per ISA tier) via
// the ikafssn_pfd_{avx2,avx512bw,avx512vbmi2} OBJECT libraries declared
// in the top-level CMakeLists.txt.  Each compilation is given:
//
//   -DFastPForLib=FastPForLib_<tier>     (renames FastPFor's namespace
//                                          at preprocessor time so the
//                                          three sets of bitpacker
//                                          symbols do not collide at
//                                          link time)
//   -DIKAFSSN_PFD_TIER_NAME=<tier>       (used here to name the
//                                          tier-specific ikafssn::pfd
//                                          inner namespace)
//   -m<arch> ...                          (controls the instructions the
//                                          bitpacker primitives and the
//                                          surrounding scalar code in
//                                          this file are allowed to use)
//
// We no longer use FastPFor's CompositeCodec / SIMDFastPFor / VByte tail.
// Instead we drive the low-level FastPForLib::simdpackwithoutmask /
// simdunpack primitives ourselves and wrap them in a custom posting
// layout that supports FOR-within-block on .kpx absolute positions.
//
// Posting layout (per k-mer, on disk, byte-aligned):
//
//   .kix v5 (delta-encoded seq_id stream):
//     [u32 count]
//     [u32 first_id]                       — only present when count >= 1
//     repeated (count - 1) / 128 times:
//       [u8 b]                             — bit-width 0..32
//       [body: 128 * b / 8 bytes]          — bitpacked deltas
//     [u8 tail_count]                       — (count - 1) % 128
//     [tail varint stream of tail_count seq_id deltas]
//
//   .kpx v5 (FOR-within-block on absolute positions):
//     [u32 count]
//     repeated count / 128 times:
//       [u8 b]                             — bit-width 0..32
//       [u32 min]                          — block FOR base
//       [body: 128 * b / 8 bytes]          — bitpacked (value - min)
//     [u8 tail_count]                       — count % 128
//     [u32 tail_min]                        — only present when tail_count >= 1
//     [tail varint stream of (value - tail_min)]
//
// The "magic" inside each tier OBJECT library is that the bitpacker
// primitives themselves are recompiled per ISA tier, so that the
// surrounding scalar loops (FOR computation, varint emit, etc.) here
// also pick up auto-vectorisation gains for AVX-512 hosts.

#include "index/pfd_codec.hpp"

// FastPForLib is rewritten by the build system to FastPForLib_<tier>;
// we only use the low-level bitpacker primitives.
#include "simdbitpacking.h"
#include "common.h"

#include <algorithm>
#include <cstring>
#include <vector>

#ifndef IKAFSSN_PFD_TIER_NAME
#error "IKAFSSN_PFD_TIER_NAME must be set per tier (avx2 / avx512bw / avx512vbmi2)"
#endif

#define IKAFSSN_PFD_TIER_NS_(x) ikafssn_pfd_##x
#define IKAFSSN_PFD_TIER_NS(x)  IKAFSSN_PFD_TIER_NS_(x)

namespace ikafssn::pfd::IKAFSSN_PFD_TIER_NS(IKAFSSN_PFD_TIER_NAME) {

namespace {

namespace pfor_ns = FastPForLib;

constexpr int kBlockSize = 128;
constexpr int kBlockAlign = 16; // __m128i

// Number of bits required to represent v.  Returns 0 for v == 0.
inline std::uint8_t bits_required(std::uint32_t v) {
    return v == 0 ? std::uint8_t(0) : std::uint8_t(32 - __builtin_clz(v));
}

// LEB128 varint encode/decode (small inline helpers, kept local to keep
// the per-tier TU self-contained).
inline std::size_t varint_encode(std::uint32_t v, std::uint8_t* dst) {
    std::size_t n = 0;
    do {
        std::uint8_t byte = v & 0x7F;
        v >>= 7;
        if (v != 0) byte |= 0x80;
        dst[n++] = byte;
    } while (v != 0);
    return n;
}

inline std::size_t varint_decode(const std::uint8_t* src, std::uint32_t& out) {
    out = 0;
    unsigned shift = 0;
    std::size_t n = 0;
    std::uint8_t byte;
    do {
        byte = src[n++];
        out |= std::uint32_t(byte & 0x7F) << shift;
        shift += 7;
    } while (byte & 0x80);
    return n;
}

// Encode one block of 128 uint32 values without FOR (used by .kix; the
// caller has already converted the seq_id stream to deltas).  Writes
// 1 + 16*b bytes to out and returns the byte count.
std::size_t encode_block_plain(const std::uint32_t* in_128,
                               std::vector<std::uint8_t>& out) {
    std::uint32_t mx = 0;
    for (int i = 0; i < kBlockSize; i++) {
        if (in_128[i] > mx) mx = in_128[i];
    }
    const std::uint8_t b = bits_required(mx);
    const std::size_t before = out.size();

    out.resize(before + 1);
    out[before] = b;
    if (b == 0) return 1;

    const std::size_t body_bytes = std::size_t(16) * b;  // 128*b/8
    out.resize(before + 1 + body_bytes);

    alignas(kBlockAlign) __m128i tmp[32];  // 32 m128i covers b up to 32
    pfor_ns::simdpackwithoutmask(in_128, tmp, b);
    std::memcpy(out.data() + before + 1, tmp, body_bytes);
    return 1 + body_bytes;
}

// Decode one block of 128 uint32 values without FOR.  Returns bytes
// consumed.
std::size_t decode_block_plain(const std::uint8_t* in,
                               std::uint32_t* out_128) {
    const std::uint8_t b = in[0];
    if (b == 0) {
        std::memset(out_128, 0, kBlockSize * sizeof(std::uint32_t));
        return 1;
    }
    const std::size_t body_bytes = std::size_t(16) * b;
    alignas(kBlockAlign) __m128i tmp[32];
    std::memcpy(tmp, in + 1, body_bytes);
    pfor_ns::simdunpack(tmp, out_128, b);
    return 1 + body_bytes;
}

// Encode one block of 128 uint32 values with FOR-within-block (used by
// .kpx).  Subtracts the block's min before bitpacking; stores min as a
// 4-byte field so the per-block bit width depends on the spread within
// the block, not on absolute value magnitude.  Writes 5 + 16*b bytes
// and returns the byte count.
std::size_t encode_block_for(const std::uint32_t* in_128,
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

    out.resize(before + 5);
    out[before] = b;
    std::memcpy(out.data() + before + 1, &mn, sizeof(std::uint32_t));
    if (b == 0) return 5;

    alignas(kBlockAlign) std::uint32_t shifted[kBlockSize];
    for (int i = 0; i < kBlockSize; i++) shifted[i] = in_128[i] - mn;

    const std::size_t body_bytes = std::size_t(16) * b;
    out.resize(before + 5 + body_bytes);

    alignas(kBlockAlign) __m128i tmp[32];
    pfor_ns::simdpackwithoutmask(shifted, tmp, b);
    std::memcpy(out.data() + before + 5, tmp, body_bytes);
    return 5 + body_bytes;
}

// Decode one FOR block; returns bytes consumed.
std::size_t decode_block_for(const std::uint8_t* in, std::uint32_t* out_128) {
    const std::uint8_t b = in[0];
    std::uint32_t mn;
    std::memcpy(&mn, in + 1, sizeof(std::uint32_t));
    if (b == 0) {
        for (int i = 0; i < kBlockSize; i++) out_128[i] = mn;
        return 5;
    }
    const std::size_t body_bytes = std::size_t(16) * b;
    alignas(kBlockAlign) __m128i tmp[32];
    std::memcpy(tmp, in + 1 + 4, body_bytes);
    pfor_ns::simdunpack(tmp, out_128, b);
    for (int i = 0; i < kBlockSize; i++) out_128[i] += mn;
    return 5 + body_bytes;
}

} // anonymous namespace

// ===== .kix encode (delta stream + plain blocks + varint tail) ======

std::size_t encode_posting_kix(const std::uint32_t* delta_array,
                               std::uint32_t count,
                               std::vector<std::uint8_t>& out) {
    const std::size_t before = out.size();
    out.resize(before + sizeof(std::uint32_t));
    std::memcpy(out.data() + before, &count, sizeof(std::uint32_t));
    if (count == 0) return out.size() - before;

    // first_id is delta_array[0] = the absolute first seq_id.
    out.resize(out.size() + sizeof(std::uint32_t));
    std::memcpy(out.data() + before + sizeof(std::uint32_t),
                &delta_array[0], sizeof(std::uint32_t));

    // The remaining (count - 1) values are deltas[1..count]; bitpack as
    // many 128-element blocks as possible, then varint the tail.
    const std::uint32_t deltas_to_pack = count - 1;
    const std::uint32_t num_blocks = deltas_to_pack / kBlockSize;
    const std::uint32_t tail_count = deltas_to_pack % kBlockSize;

    for (std::uint32_t b = 0; b < num_blocks; b++) {
        encode_block_plain(delta_array + 1 + b * kBlockSize, out);
    }

    // Tail header + varint stream.
    out.push_back(static_cast<std::uint8_t>(tail_count));
    if (tail_count > 0) {
        std::uint8_t tmp[5];
        for (std::uint32_t i = 0; i < tail_count; i++) {
            const std::uint32_t v = delta_array[1 + num_blocks * kBlockSize + i];
            const std::size_t n = varint_encode(v, tmp);
            out.insert(out.end(), tmp, tmp + n);
        }
    }
    return out.size() - before;
}

// ===== .kpx encode (FOR blocks + varint tail with FOR base) =========

std::size_t encode_posting_kpx(const std::uint32_t* abs_pos_array,
                               std::uint32_t count,
                               std::vector<std::uint8_t>& out) {
    const std::size_t before = out.size();
    out.resize(before + sizeof(std::uint32_t));
    std::memcpy(out.data() + before, &count, sizeof(std::uint32_t));
    if (count == 0) return out.size() - before;

    const std::uint32_t num_blocks = count / kBlockSize;
    const std::uint32_t tail_count = count % kBlockSize;

    for (std::uint32_t b = 0; b < num_blocks; b++) {
        encode_block_for(abs_pos_array + b * kBlockSize, out);
    }

    out.push_back(static_cast<std::uint8_t>(tail_count));
    if (tail_count > 0) {
        // FOR-encode the tail too: emit min + varint of (value - min).
        const std::uint32_t* tail = abs_pos_array + num_blocks * kBlockSize;
        std::uint32_t mn = tail[0];
        for (std::uint32_t i = 1; i < tail_count; i++) {
            if (tail[i] < mn) mn = tail[i];
        }
        out.resize(out.size() + sizeof(std::uint32_t));
        std::memcpy(out.data() + out.size() - sizeof(std::uint32_t),
                    &mn, sizeof(std::uint32_t));

        std::uint8_t tmp[5];
        for (std::uint32_t i = 0; i < tail_count; i++) {
            const std::size_t n = varint_encode(tail[i] - mn, tmp);
            out.insert(out.end(), tmp, tmp + n);
        }
    }
    return out.size() - before;
}

// ===== open_stream: decode the entire posting into the StreamCtx ====

bool open_stream_kix(const std::uint8_t* posting, std::size_t bytes,
                     ikafssn::pfd::StreamCtx& ctx) {
    ctx.decoded.clear();
    ctx.count = 0;
    ctx.pos = 0;
    if (bytes == 0) return true;
    if (bytes < sizeof(std::uint32_t)) return false;

    std::uint32_t count;
    std::memcpy(&count, posting, sizeof(std::uint32_t));
    if (count == 0) return true;
    if (bytes < 2 * sizeof(std::uint32_t)) return false;

    std::uint32_t first_id;
    std::memcpy(&first_id, posting + sizeof(std::uint32_t),
                sizeof(std::uint32_t));

    ctx.decoded.resize(count);
    ctx.decoded[0] = first_id;

    const std::uint32_t deltas_to_unpack = count - 1;
    const std::uint32_t num_blocks = deltas_to_unpack / kBlockSize;
    const std::uint32_t tail_count = deltas_to_unpack % kBlockSize;

    const std::uint8_t* p = posting + 2 * sizeof(std::uint32_t);
    const std::uint8_t* end = posting + bytes;

    // Decode each full block and accumulate deltas back to absolute IDs.
    for (std::uint32_t b = 0; b < num_blocks; b++) {
        if (p >= end) return false;
        alignas(kBlockAlign) std::uint32_t deltas[kBlockSize];
        const std::size_t n = decode_block_plain(p, deltas);
        if (p + n > end) return false;
        p += n;
        std::uint32_t prev = ctx.decoded[1 + b * kBlockSize - 1];
        std::uint32_t* out = ctx.decoded.data() + 1 + b * kBlockSize;
        for (int i = 0; i < kBlockSize; i++) {
            prev += deltas[i];
            out[i] = prev;
        }
    }

    if (p >= end) return false;
    const std::uint8_t got_tail = *p++;
    if (got_tail != tail_count) return false;

    if (tail_count > 0) {
        std::uint32_t prev = ctx.decoded[count - 1 - tail_count];
        for (std::uint32_t i = 0; i < tail_count; i++) {
            if (p >= end) return false;
            std::uint32_t d;
            const std::size_t n = varint_decode(p, d);
            if (p + n > end) return false;
            p += n;
            prev += d;
            ctx.decoded[count - tail_count + i] = prev;
        }
    }

    ctx.count = count;
    ctx.pos = 0;
    return true;
}

bool open_stream_kpx(const std::uint8_t* posting, std::size_t bytes,
                     ikafssn::pfd::StreamCtx& ctx) {
    ctx.decoded.clear();
    ctx.count = 0;
    ctx.pos = 0;
    if (bytes == 0) return true;
    if (bytes < sizeof(std::uint32_t)) return false;

    std::uint32_t count;
    std::memcpy(&count, posting, sizeof(std::uint32_t));
    if (count == 0) return true;

    ctx.decoded.resize(count);

    const std::uint32_t num_blocks = count / kBlockSize;
    const std::uint32_t tail_count = count % kBlockSize;

    const std::uint8_t* p = posting + sizeof(std::uint32_t);
    const std::uint8_t* end = posting + bytes;

    for (std::uint32_t b = 0; b < num_blocks; b++) {
        if (p >= end) return false;
        const std::size_t n = decode_block_for(p, ctx.decoded.data() + b * kBlockSize);
        if (p + n > end) return false;
        p += n;
    }

    if (p >= end) return false;
    const std::uint8_t got_tail = *p++;
    if (got_tail != tail_count) return false;

    if (tail_count > 0) {
        if (p + sizeof(std::uint32_t) > end) return false;
        std::uint32_t tail_min;
        std::memcpy(&tail_min, p, sizeof(std::uint32_t));
        p += sizeof(std::uint32_t);

        std::uint32_t* tail_out = ctx.decoded.data() + num_blocks * kBlockSize;
        for (std::uint32_t i = 0; i < tail_count; i++) {
            if (p >= end) return false;
            std::uint32_t d;
            const std::size_t n = varint_decode(p, d);
            if (p + n > end) return false;
            p += n;
            tail_out[i] = tail_min + d;
        }
    }

    ctx.count = count;
    ctx.pos = 0;
    return true;
}

} // namespace ikafssn::pfd::ikafssn_pfd_<tier>
