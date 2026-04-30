// Phase 5i — per-tier custom block codec.
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
// For .kix postings (sorted distinct seq_id delta stream — v7) we drive
// FastPFor's CompositeCodec<SIMDFastPFor<4>, VariableByte> directly.
//
// For .kpx postings (absolute position stream) Phase 5g-2 splits each
// posting into per-(kmer, seq_id) partition groups (when the within-
// sequence occurrence count exceeds a build-time threshold) plus one
// short bucket for the remaining low-multiplicity occurrences.  Phase 5i
// makes the short bucket self-describing (its own delta-encoded seq_id
// list + per-seq_id u8 occurrence counts) so decoding is candidate-set-
// driven and independent of the .kix stream.

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

inline std::uint8_t bits_required(std::uint32_t v) {
    return v == 0 ? std::uint8_t(0) : std::uint8_t(32 - __builtin_clz(v));
}

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

// ===== .kix v7 encode (distinct seq_id delta stream → SIMDFastPFor + VByte tail) =====
//
// On-disk per-posting blob:
//   [u32 distinct_count]  — number of distinct seq_ids represented
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

// ===== .kpx v7 encode (per-(kmer, seq_id) partitioning + self-describing short bucket) =====

std::size_t encode_posting_kpx(const std::uint32_t* distinct_sid,
                               const std::uint8_t*  occ_count,
                               std::uint32_t distinct_count,
                               const std::uint32_t* abs_pos_array,
                               std::uint32_t position_count,
                               std::uint32_t freq_threshold_part,
                               std::vector<std::uint8_t>& out) {
    const std::size_t before = out.size();

    // Header placeholders: distinct_count, position_count, partition_count.
    out.resize(before + 3 * sizeof(std::uint32_t));
    std::memcpy(out.data() + before, &distinct_count, sizeof(std::uint32_t));
    std::memcpy(out.data() + before + sizeof(std::uint32_t),
                &position_count, sizeof(std::uint32_t));

    if (distinct_count == 0) {
        std::uint32_t zero = 0;
        std::memcpy(out.data() + before + 2 * sizeof(std::uint32_t),
                    &zero, sizeof(std::uint32_t));
        // short_seq_count = 0, short_position_count = 0
        out.resize(out.size() + 2 * sizeof(std::uint32_t));
        std::memcpy(out.data() + out.size() - 2 * sizeof(std::uint32_t),
                    &zero, sizeof(std::uint32_t));
        std::memcpy(out.data() + out.size() - sizeof(std::uint32_t),
                    &zero, sizeof(std::uint32_t));
        return out.size() - before;
    }

    // Pass 1: classify each distinct seq_id into partition vs short bucket
    // based on occ_count[k] > freq_threshold_part.
    std::uint32_t partition_count = 0;
    std::uint32_t short_seq_count = 0;
    std::uint32_t short_position_count = 0;
    for (std::uint32_t k = 0; k < distinct_count; k++) {
        if (occ_count[k] > freq_threshold_part) {
            partition_count++;
        } else {
            short_seq_count++;
            short_position_count += occ_count[k];
        }
    }

    std::memcpy(out.data() + before + 2 * sizeof(std::uint32_t),
                &partition_count, sizeof(std::uint32_t));

    // Pass 2: emit partition groups in seq_id order.
    std::uint32_t pos_cursor = 0;
    for (std::uint32_t k = 0; k < distinct_count; k++) {
        const std::uint32_t cnt = occ_count[k];
        if (occ_count[k] > freq_threshold_part) {
            const std::size_t hdr_off = out.size();
            out.resize(hdr_off + 2 * sizeof(std::uint32_t));
            std::memcpy(out.data() + hdr_off,
                        &distinct_sid[k], sizeof(std::uint32_t));
            std::memcpy(out.data() + hdr_off + sizeof(std::uint32_t),
                        &cnt, sizeof(std::uint32_t));
            encode_for_stream(abs_pos_array + pos_cursor, cnt, out);
        }
        pos_cursor += cnt;
    }

    // Pass 3: emit short bucket header + delta seq_ids + occ_counts + positions.
    {
        const std::size_t cnt_off = out.size();
        out.resize(cnt_off + 2 * sizeof(std::uint32_t));
        std::memcpy(out.data() + cnt_off,
                    &short_seq_count, sizeof(std::uint32_t));
        std::memcpy(out.data() + cnt_off + sizeof(std::uint32_t),
                    &short_position_count, sizeof(std::uint32_t));
    }

    if (short_seq_count > 0) {
        // Delta-encoded seq_id list: first absolute, then varint deltas.
        bool first = true;
        std::uint32_t prev_sid = 0;
        std::uint8_t tmp[5];
        for (std::uint32_t k = 0; k < distinct_count; k++) {
            if (occ_count[k] > freq_threshold_part) continue;
            const std::uint32_t sid = distinct_sid[k];
            if (first) {
                out.resize(out.size() + sizeof(std::uint32_t));
                std::memcpy(out.data() + out.size() - sizeof(std::uint32_t),
                            &sid, sizeof(std::uint32_t));
                first = false;
            } else {
                const std::uint32_t d = sid - prev_sid;
                const std::size_t n = varint_encode(d, tmp);
                out.insert(out.end(), tmp, tmp + n);
            }
            prev_sid = sid;
        }

        // Per-seq_id u8 occ_count list (in the same order as the seq_id list).
        for (std::uint32_t k = 0; k < distinct_count; k++) {
            if (occ_count[k] > freq_threshold_part) continue;
            out.push_back(occ_count[k]);
        }

        // Concatenated positions (intra-seq ordering preserved).
        std::vector<std::uint32_t> short_buf;
        short_buf.reserve(short_position_count);
        std::uint32_t cursor2 = 0;
        for (std::uint32_t k = 0; k < distinct_count; k++) {
            const std::uint32_t cnt = occ_count[k];
            if (occ_count[k] <= freq_threshold_part) {
                short_buf.insert(short_buf.end(),
                                 abs_pos_array + cursor2,
                                 abs_pos_array + cursor2 + cnt);
            }
            cursor2 += cnt;
        }
        encode_for_stream(short_buf.data(), short_position_count, out);
    }

    return out.size() - before;
}

// ===== open_stream_kix: decode the entire .kix posting into the StreamCtx =====

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

    std::vector<std::uint32_t> codec_in(payload_words);
    std::memcpy(codec_in.data(), posting + 8, payload_bytes);

    ctx.decoded.resize(count);
    std::size_t nvalue = ctx.decoded.size();
    kix_codec().decodeArray(codec_in.data(), payload_words,
                            ctx.decoded.data(), nvalue);
    if (nvalue != count) return false;

    // v7: encoder writes [abs_first, d1, d2, ...] over the **distinct**
    // seq_id stream; cumulative sum reconstructs absolute distinct seq_ids.
    for (std::uint32_t i = 1; i < count; i++) {
        ctx.decoded[i] += ctx.decoded[i - 1];
    }

    ctx.count = count;
    ctx.pos = 0;
    return true;
}

// ===== open_stream_kpx_for_candidates: candidate-set-driven decode =====

bool open_stream_kpx_for_candidates(
        const std::uint8_t* posting, std::size_t bytes,
        const std::uint32_t* candidates, std::size_t n_candidates,
        std::vector<std::vector<std::uint32_t>>& out) {

    out.assign(n_candidates, std::vector<std::uint32_t>{});

    if (n_candidates == 0) return true;
    if (bytes == 0) return true;
    if (bytes < 5 * sizeof(std::uint32_t)) return false;

    const std::uint8_t* p = posting;
    const std::uint8_t* end = posting + bytes;

    std::uint32_t distinct_count;
    std::uint32_t position_count;
    std::uint32_t partition_count;
    std::memcpy(&distinct_count,  p, sizeof(std::uint32_t)); p += sizeof(std::uint32_t);
    std::memcpy(&position_count,  p, sizeof(std::uint32_t)); p += sizeof(std::uint32_t);
    std::memcpy(&partition_count, p, sizeof(std::uint32_t)); p += sizeof(std::uint32_t);
    (void)position_count;

    if (distinct_count == 0) {
        // Expect short_seq_count = 0 + short_position_count = 0 trailer.
        if (p + 2 * sizeof(std::uint32_t) > end) return false;
        std::uint32_t a, b;
        std::memcpy(&a, p, sizeof(std::uint32_t)); p += sizeof(std::uint32_t);
        std::memcpy(&b, p, sizeof(std::uint32_t)); p += sizeof(std::uint32_t);
        if (a != 0 || b != 0) return false;
        return true;
    }

    // Helper: locate candidate index for seq_id sid via binary search.
    // Returns size_t(-1) when sid is not in the candidate list.  The two
    // streams (partition groups and short bucket) are each sorted, but
    // their union interleaves arbitrarily against the candidate list, so
    // a shared monotonic cursor is not sufficient — binary search keeps
    // both passes correct without extra bookkeeping.
    auto find_candidate = [&](std::uint32_t sid) -> std::size_t {
        std::size_t lo = 0;
        std::size_t hi = n_candidates;
        while (lo < hi) {
            std::size_t mid = lo + (hi - lo) / 2;
            if (candidates[mid] < sid)      lo = mid + 1;
            else if (candidates[mid] > sid) hi = mid;
            else return mid;
        }
        return static_cast<std::size_t>(-1);
    };

    // Pass 1: partition groups.
    for (std::uint32_t g = 0; g < partition_count; g++) {
        if (p + 2 * sizeof(std::uint32_t) > end) return false;
        std::uint32_t gsid, gcnt;
        std::memcpy(&gsid, p, sizeof(std::uint32_t)); p += sizeof(std::uint32_t);
        std::memcpy(&gcnt, p, sizeof(std::uint32_t)); p += sizeof(std::uint32_t);
        if (gcnt == 0) return false;

        const std::size_t cidx = find_candidate(gsid);
        if (cidx != static_cast<std::size_t>(-1)) {
            std::vector<std::uint32_t> buf(gcnt);
            if (!decode_for_stream(p, end, gcnt, buf.data())) return false;
            out[cidx] = std::move(buf);
        } else {
            std::vector<std::uint32_t> scratch(gcnt);
            if (!decode_for_stream(p, end, gcnt, scratch.data())) return false;
        }
    }

    // Pass 2: short bucket.
    if (p + 2 * sizeof(std::uint32_t) > end) return false;
    std::uint32_t short_seq_count;
    std::uint32_t short_position_count;
    std::memcpy(&short_seq_count,      p, sizeof(std::uint32_t)); p += sizeof(std::uint32_t);
    std::memcpy(&short_position_count, p, sizeof(std::uint32_t)); p += sizeof(std::uint32_t);

    if (short_seq_count == 0) {
        if (short_position_count != 0) return false;
        return true;
    }

    if (p + sizeof(std::uint32_t) > end) return false;
    std::vector<std::uint32_t> short_sid(short_seq_count);
    std::memcpy(&short_sid[0], p, sizeof(std::uint32_t));
    p += sizeof(std::uint32_t);
    for (std::uint32_t i = 1; i < short_seq_count; i++) {
        if (p >= end) return false;
        std::uint32_t d;
        const std::size_t n = varint_decode(p, d);
        if (p + n > end) return false;
        p += n;
        short_sid[i] = short_sid[i - 1] + d;
    }

    if (p + short_seq_count > end) return false;
    const std::uint8_t* short_occ = p;
    p += short_seq_count;

    std::vector<std::uint32_t> short_positions(short_position_count);
    if (short_position_count > 0) {
        if (!decode_for_stream(p, end, short_position_count, short_positions.data())) {
            return false;
        }
    }

    std::size_t pos_cursor = 0;
    for (std::uint32_t k = 0; k < short_seq_count; k++) {
        const std::uint32_t sid = short_sid[k];
        const std::uint32_t cnt = short_occ[k];
        const std::size_t cidx = find_candidate(sid);
        if (cidx != static_cast<std::size_t>(-1)) {
            out[cidx].assign(short_positions.begin() + pos_cursor,
                             short_positions.begin() + pos_cursor + cnt);
        }
        pos_cursor += cnt;
    }

    if (pos_cursor != short_position_count) return false;
    return true;
}

} // namespace ikafssn::pfd::ikafssn_pfd_<tier>
