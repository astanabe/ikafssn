// Phase 5g-2 — per-tier custom block codec.
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
// For .kpx postings (absolute position stream) Phase 5g-2 splits each
// posting into per-(kmer, seq_id) partition groups (when the within-
// sequence occurrence count exceeds a build-time threshold) plus one
// short bucket for the remaining low-multiplicity occurrences.  Each
// stream uses the Phase 5e FOR-within-block layout — encode_block_for /
// decode_block_for here are unchanged from v5.

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

// Encode `count` absolute positions through the FOR-block stream
// (count/128 full blocks followed by an optional tail).  Identical to the
// Phase 5e payload that encode_posting_kpx used to emit, factored out
// here so v6 can reuse it for both partition groups and the short bucket.
void encode_for_stream(const std::uint32_t* abs_pos_array,
                       std::uint32_t count,
                       std::vector<std::uint8_t>& out) {
    const std::uint32_t num_blocks = count / kBlockSize;
    const std::uint32_t tail_count = count % kBlockSize;

    for (std::uint32_t b = 0; b < num_blocks; b++) {
        encode_block_for(abs_pos_array + b * kBlockSize, out);
    }

    out.push_back(static_cast<std::uint8_t>(tail_count));
    if (tail_count > 0) {
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
}

// Decode a FOR-block stream of `count` positions into `out`.  Returns
// false if the byte cursor `p` would run past `end` (corrupt index);
// `p` is advanced to the first byte after the stream on success.
bool decode_for_stream(const std::uint8_t*& p, const std::uint8_t* end,
                       std::uint32_t count, std::uint32_t* out) {
    const std::uint32_t num_blocks = count / kBlockSize;
    const std::uint32_t tail_count = count % kBlockSize;

    for (std::uint32_t b = 0; b < num_blocks; b++) {
        if (p >= end) return false;
        const std::size_t n = decode_block_for(p, out + b * kBlockSize);
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

        std::uint32_t* tail_out = out + num_blocks * kBlockSize;
        for (std::uint32_t i = 0; i < tail_count; i++) {
            if (p >= end) return false;
            std::uint32_t d;
            const std::size_t n = varint_decode(p, d);
            if (p + n > end) return false;
            p += n;
            tail_out[i] = tail_min + d;
        }
    }
    return true;
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

// ===== .kpx encode (v6: per-(kmer, seq_id) partitioning + short bucket) =====

std::size_t encode_posting_kpx(const std::uint32_t* sid_array,
                               const std::uint32_t* abs_pos_array,
                               std::uint32_t count,
                               std::uint32_t freq_threshold_part,
                               std::vector<std::uint8_t>& out) {
    const std::size_t before = out.size();

    // Header placeholders: total_count, partition_count.
    out.resize(before + 2 * sizeof(std::uint32_t));
    std::memcpy(out.data() + before, &count, sizeof(std::uint32_t));

    if (count == 0) {
        std::uint32_t zero = 0;
        std::memcpy(out.data() + before + sizeof(std::uint32_t),
                    &zero, sizeof(std::uint32_t));
        // short_bucket_count = 0
        out.resize(out.size() + sizeof(std::uint32_t));
        std::memcpy(out.data() + out.size() - sizeof(std::uint32_t),
                    &zero, sizeof(std::uint32_t));
        return out.size() - before;
    }

    // Pass 1: identify partition-group runs in sid_array (sorted, may
    // contain duplicates).  A run with length > freq_threshold_part
    // becomes a partition group; the rest go into the short bucket.
    struct RunInfo {
        std::uint32_t seq_id;
        std::uint32_t start;   // index into sid_array / abs_pos_array
        std::uint32_t length;
    };
    std::vector<RunInfo> runs;
    runs.reserve(64);
    {
        std::uint32_t i = 0;
        while (i < count) {
            std::uint32_t j = i + 1;
            while (j < count && sid_array[j] == sid_array[i]) j++;
            runs.push_back({sid_array[i], i, j - i});
            i = j;
        }
    }

    std::uint32_t partition_count = 0;
    std::uint32_t short_count = 0;
    for (const auto& r : runs) {
        if (r.length > freq_threshold_part) {
            partition_count++;
        } else {
            short_count += r.length;
        }
    }

    std::memcpy(out.data() + before + sizeof(std::uint32_t),
                &partition_count, sizeof(std::uint32_t));

    // Pass 2: emit partition groups in seq_id order (runs are already
    // ordered because sid_array is sorted).
    for (const auto& r : runs) {
        if (r.length <= freq_threshold_part) continue;

        const std::size_t hdr_off = out.size();
        out.resize(hdr_off + 2 * sizeof(std::uint32_t));
        std::memcpy(out.data() + hdr_off,
                    &r.seq_id, sizeof(std::uint32_t));
        std::memcpy(out.data() + hdr_off + sizeof(std::uint32_t),
                    &r.length, sizeof(std::uint32_t));
        encode_for_stream(abs_pos_array + r.start, r.length, out);
    }

    // Pass 3: emit short bucket — concat all sub-threshold runs in input
    // order (which is the seq_id-sorted order, identical to the order the
    // builder is going to encounter them at decode time when scanning
    // sid_stream).
    {
        const std::size_t cnt_off = out.size();
        out.resize(cnt_off + sizeof(std::uint32_t));
        std::memcpy(out.data() + cnt_off, &short_count, sizeof(std::uint32_t));
    }

    if (short_count > 0) {
        std::vector<std::uint32_t> short_buf;
        short_buf.reserve(short_count);
        for (const auto& r : runs) {
            if (r.length > freq_threshold_part) continue;
            short_buf.insert(short_buf.end(),
                             abs_pos_array + r.start,
                             abs_pos_array + r.start + r.length);
        }
        encode_for_stream(short_buf.data(), short_count, out);
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
                     const std::uint32_t* sid_stream, std::size_t n_sids,
                     ikafssn::pfd::StreamCtx& ctx) {
    ctx.decoded.clear();
    ctx.count = 0;
    ctx.pos = 0;
    if (bytes == 0) return true;
    if (bytes < 3 * sizeof(std::uint32_t)) return false;

    const std::uint8_t* p = posting;
    const std::uint8_t* end = posting + bytes;

    std::uint32_t total_count;
    std::uint32_t partition_count;
    std::memcpy(&total_count,     p, sizeof(std::uint32_t));
    p += sizeof(std::uint32_t);
    std::memcpy(&partition_count, p, sizeof(std::uint32_t));
    p += sizeof(std::uint32_t);

    if (total_count == 0) {
        // Expect short_bucket_count = 0 trailer.
        if (p + sizeof(std::uint32_t) > end) return false;
        std::uint32_t sbc;
        std::memcpy(&sbc, p, sizeof(std::uint32_t));
        if (sbc != 0) return false;
        return true;
    }

    if (sid_stream == nullptr || n_sids != total_count) return false;

    // Decode each partition group into a temporary buffer; remember the
    // per-group seq_id and length so the lock-step merge below can index
    // into it without re-scanning the byte stream.
    struct GroupView {
        std::uint32_t seq_id;
        std::uint32_t count;
        std::uint32_t buf_start; // index into all_partition_buf
    };
    std::vector<GroupView> groups;
    groups.reserve(partition_count);
    std::vector<std::uint32_t> all_partition_buf;

    for (std::uint32_t g = 0; g < partition_count; g++) {
        if (p + 2 * sizeof(std::uint32_t) > end) return false;
        std::uint32_t gsid, gcnt;
        std::memcpy(&gsid, p, sizeof(std::uint32_t));
        p += sizeof(std::uint32_t);
        std::memcpy(&gcnt, p, sizeof(std::uint32_t));
        p += sizeof(std::uint32_t);
        if (gcnt == 0) return false; // partition groups are non-empty

        const std::size_t buf_start = all_partition_buf.size();
        all_partition_buf.resize(buf_start + gcnt);
        if (!decode_for_stream(p, end, gcnt, all_partition_buf.data() + buf_start)) {
            return false;
        }
        groups.push_back({gsid, gcnt, static_cast<std::uint32_t>(buf_start)});
    }

    if (p + sizeof(std::uint32_t) > end) return false;
    std::uint32_t short_count;
    std::memcpy(&short_count, p, sizeof(std::uint32_t));
    p += sizeof(std::uint32_t);

    std::vector<std::uint32_t> short_buf(short_count);
    if (short_count > 0) {
        if (!decode_for_stream(p, end, short_count, short_buf.data())) {
            return false;
        }
    }

    // Merge into ctx.decoded[i] aligned with sid_stream[i].
    ctx.decoded.resize(total_count);

    std::uint32_t pi = 0;
    std::uint32_t group_cursor = 0;
    std::uint32_t short_cursor = 0;
    for (std::uint32_t k = 0; k < total_count; k++) {
        const std::uint32_t sid = sid_stream[k];
        if (pi < partition_count && groups[pi].seq_id == sid) {
            ctx.decoded[k] = all_partition_buf[groups[pi].buf_start + group_cursor];
            group_cursor++;
            if (group_cursor == groups[pi].count) {
                pi++;
                group_cursor = 0;
            }
        } else {
            if (short_cursor >= short_count) return false;
            ctx.decoded[k] = short_buf[short_cursor++];
        }
    }

    if (pi != partition_count) return false;
    if (short_cursor != short_count) return false;

    ctx.count = total_count;
    ctx.pos = 0;
    return true;
}

} // namespace ikafssn::pfd::ikafssn_pfd_<tier>
