// Phase 5g-1 — per-tier custom block codec.
//
// This translation unit is compiled FOUR TIMES (once per ISA tier) via
// the ikafssn_pfd_{sse42,avx2,avx512bw,avx512vbmi2} OBJECT libraries
// declared in the top-level CMakeLists.txt.  Each compilation is given:
//
//   -DFastPForLib=FastPForLib_<tier>     (renames FastPFor's namespace
//                                          at preprocessor time so the
//                                          four sets of symbols do not
//                                          collide at link time)
//   -DIKAFSSN_PFD_TIER_NAME=<tier>       (used here to name the
//                                          tier-specific ikafssn::pfd
//                                          inner namespace)
//   -m<arch> ...                          (controls the instructions the
//                                          bitpacker primitives and the
//                                          surrounding scalar code in
//                                          this file are allowed to use)
//
// For .kix postings (sorted seq_id delta stream) we drive FastPFor's
// CompositeCodec<SIMDFastPFor<4>, VariableByte> directly — that's the
// classical PForDelta layout (per-block bit-width with an out-of-line
// VByte exception stream for the few values that exceed it), which gives
// us back the patched-PFor advantage that Phase 5e accidentally lost.
//
// For .kpx postings (absolute position stream) we keep the Phase 5e
// custom FOR-within-block layout: each 128-element block subtracts its
// min before bitpacking so the per-block bit-width depends on the
// within-sequence spread rather than absolute position magnitude.
// FastPFor's CompositeCodec doesn't support per-block FOR base, so we
// continue to drive simdpackwithoutmask / simdunpack directly here.

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

// FastPFor codec used for .kix postings (PForDelta + VByte exception
// stream).  Stateless aside from a small per-codec scratch buffer, but
// the encoder/decoder methods are not const, so we keep a thread_local
// instance to avoid contention without paying construction overhead on
// every posting.
using KixCodec = pfor_ns::CompositeCodec<
    pfor_ns::SIMDFastPFor<4>,
    pfor_ns::VariableByte>;

KixCodec& kix_codec() {
    thread_local KixCodec codec;
    return codec;
}

// LEB128 varint helpers (used by the .kpx tail stream).
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

// Number of bits required to represent v.  Returns 0 for v == 0.
inline std::uint8_t bits_required(std::uint32_t v) {
    return v == 0 ? std::uint8_t(0) : std::uint8_t(32 - __builtin_clz(v));
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

// ===== .kix encode (delta stream → SIMDFastPFor + VByte tail) =======
//
// The builder passes `[abs_first, d1, d2, ..., d_{count-1}]`, i.e. the
// first absolute seq_id followed by consecutive deltas.  SIMDFastPFor's
// PForDelta layout handles the abs_first wide-value case via its
// VByte exception stream, so we feed the array straight in.  open_stream
// reconstructs absolutes by cumulative sum at decode time.
//
// On-disk per-posting blob:
//   [u32 count]           — number of seq_id values represented
//   [u32 payload_words]   — number of u32 words written by the codec
//   [u32 payload[N]]      — codec output, byte-unaligned

std::size_t encode_posting_kix(const std::uint32_t* delta_array,
                               std::uint32_t count,
                               std::vector<std::uint8_t>& out) {
    const std::size_t before = out.size();
    out.resize(before + 8);
    std::memcpy(out.data() + before, &count, sizeof(std::uint32_t));
    if (count == 0) {
        std::uint32_t zero = 0;
        std::memcpy(out.data() + before + 4, &zero, sizeof(std::uint32_t));
        return 8;
    }

    // FastPFor's encode interface is "give me a uint32_t* output buffer
    // big enough for the worst case".  Worst case = original size plus a
    // small constant for the VByte tail.  Reserve len*2 + 1024 to be
    // safe; the codec will report the actual nvalue used.
    const std::size_t worst_words = std::size_t(count) * 2 + 1024;
    std::vector<std::uint32_t> codec_out(worst_words);
    std::size_t nvalue = codec_out.size();
    kix_codec().encodeArray(delta_array, count, codec_out.data(), nvalue);

    const std::uint32_t payload_words = static_cast<std::uint32_t>(nvalue);
    std::memcpy(out.data() + before + 4, &payload_words, sizeof(std::uint32_t));

    const std::size_t payload_bytes = std::size_t(payload_words) * sizeof(std::uint32_t);
    out.resize(out.size() + payload_bytes);
    std::memcpy(out.data() + before + 8, codec_out.data(), payload_bytes);

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
    if (bytes < 8) return false;

    std::uint32_t count;
    std::uint32_t payload_words;
    std::memcpy(&count,         posting,     sizeof(std::uint32_t));
    std::memcpy(&payload_words, posting + 4, sizeof(std::uint32_t));
    if (count == 0) return true;

    const std::size_t payload_bytes = std::size_t(payload_words) * sizeof(std::uint32_t);
    if (bytes < 8 + payload_bytes) return false;

    // FastPFor's decoder reads uint32_t-aligned input; the on-wire
    // payload is byte-aligned (offsets[i] is byte-granular), so copy it
    // into an aligned scratch buffer first.
    std::vector<std::uint32_t> codec_in(payload_words);
    std::memcpy(codec_in.data(), posting + 8, payload_bytes);

    ctx.decoded.resize(count);
    std::size_t nvalue = ctx.decoded.size();
    kix_codec().decodeArray(codec_in.data(), payload_words,
                            ctx.decoded.data(), nvalue);
    if (nvalue != count) return false;

    // Builder writes the array as [abs_first, d1, d2, ...]; reconstruct
    // absolute seq_ids via cumulative sum so SeqIdDecoder consumers see
    // absolutes (the field name `delta_array` and the cumsum convention
    // are preserved from the v3/v4 layout).
    for (std::uint32_t i = 1; i < count; i++) {
        ctx.decoded[i] += ctx.decoded[i - 1];
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
