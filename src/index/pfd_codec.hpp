#pragma once

// Phase 5i — split codec for .kix v7 and .kpx v7.
//
//   .kix v7: FastPFor CompositeCodec<SIMDFastPFor<4>, VariableByte>
//            (PForDelta + VByte tail) over the **distinct seq_id**
//            delta stream.  Intra-sequence k-mer duplicates are
//            removed by a SIMD dedup kernel (src/index/seq_id_dedup.*)
//            at build time, so the input stream to the codec is
//            [abs_first, d1, d2, ...] with d_i >= 1 strictly.
//            Posting layout on disk:
//              [u32 distinct_count]         — distinct seq_ids in this k-mer
//              [u32 payload_words]          — payload size in u32 words
//              [u32 payload[payload_words]] — codec output, byte-unaligned
//                                             on the wire (use memcpy)
//
//   .kpx v7: per-(kmer, seq_id) partitioned position posting with a
//            **self-describing short bucket**.  A (k-mer, seq_id) cluster
//            whose occurrence count exceeds the build-time
//            `freq_threshold_part` (max 255) is split out into its own
//            partition group; remaining occurrences are merged into a
//            short bucket carrying its own delta-encoded seq_id list and
//            per-seq_id u8 occurrence counts.  The decoder no longer
//            needs to walk lock-step against the .kix stream — callers
//            pass a sorted candidate seq_id set instead.
//            Posting layout on disk:
//              [u32 distinct_count]         — distinct seq_ids (matches .kix)
//              [u32 position_count]         — total position count
//              [u32 partition_count]
//              repeated partition_count times (sorted by seq_id):
//                [u32 seq_id]
//                [u32 occurrence_count]
//                [FOR-block stream over occurrence_count positions]
//              [u32 short_seq_count]        — distinct seq_ids in short bucket
//              [u32 short_position_count]   — total positions in short bucket
//              if short_seq_count > 0:
//                [u32 abs_first_seq_id]
//                [varint  delta_seq_id[short_seq_count - 1]]
//                [u8      occ_count[short_seq_count]]
//                [FOR-block stream over short_position_count positions]
//            Each FOR-block stream is the same encode_block_for + tail
//            FOR layout used since Phase 5e:
//                repeated count/128 times:
//                  [u8 b][u32 min][128*b/8 bytes bitpacked (value-min)]
//                [u8 tail_count]
//                if tail_count > 0:
//                  [u32 tail_min][varint stream of (value - tail_min)]
//
// .kix open_stream materialises the entire decoded posting (absolute
// distinct seq_ids) into StreamCtx::decoded.  .kpx decoding is now
// candidate-set-driven (open_stream_kpx_for_candidates): the caller
// hands a sorted candidate seq_id list and receives per-candidate
// position vectors.

#include <cstdint>
#include <cstddef>
#include <cstring>
#include <memory>
#include <vector>

namespace ikafssn::pfd {

// FastPFor SIMD block size (the codec we wrap is simdfastpfor128 =
// CompositeCodec<SIMDFastPFor<4>, VariableByte>; per-block element count is
// 128, matching the plan).
inline constexpr int kPfdBlockSize = 128;

// Per-posting byte-stream header for .kix (distinct_count + payload_words).
// The .kpx layout uses its own variable-length header (see the file-level
// comment above) so this constant only describes the .kix posting blob.
inline constexpr size_t kPostingHeaderBytes = 8;

// === posting-level encode wrappers ===

// Encode the distinct seq_id-delta stream for a .kix posting (v7).  Writes
// distinct_count + payload_words + payload into `out` (appended).  Returns
// the number of bytes written.
size_t encode_posting_kix(const uint32_t* delta_array, uint32_t count,
                          std::vector<uint8_t>& out);

// Encode the absolute-position stream for a .kpx posting (v7).
//   distinct_sid          length = distinct_count, sorted ascending, no dups
//   occ_count             length = distinct_count, occurrence per distinct_sid
//                         (must be <= freq_threshold_part for short-bucket
//                         entries; long entries are spilled into partition
//                         groups based on the threshold)
//   distinct_count        number of distinct seq_ids
//   abs_pos_array         length = position_count; positions are grouped by
//                         seq_id (matching distinct_sid order) and sorted
//                         ascending within each group
//   position_count        sum of occ_count
//   freq_threshold_part   max 255; per-(kmer, seq_id) occurrence count
//                         strictly greater becomes a partition group
size_t encode_posting_kpx(const uint32_t* distinct_sid,
                          const uint8_t*  occ_count,
                          uint32_t distinct_count,
                          const uint32_t* abs_pos_array,
                          uint32_t position_count,
                          uint32_t freq_threshold_part,
                          std::vector<uint8_t>& out);

// === streaming decode context (for .kix only since v7) ===
//
// open_stream_kix materialises the entire decoded posting into `decoded`.
// .kpx decoding moved to a candidate-set-driven API (see below).
struct StreamCtx {
    std::vector<uint32_t> decoded;
    uint32_t count = 0;     // total elements in `decoded`
    uint32_t pos   = 0;     // next read index
};

// Initialise a StreamCtx for a .kix posting starting at `posting`.
//   posting:    pointer to first byte of the per-kmer posting blob
//   bytes:      total byte length of the posting blob
//
// Returns false on header / payload size mismatch (corrupt index).
bool open_stream_kix(const uint8_t* posting, size_t bytes, StreamCtx& ctx);

// Read up to `max_count` decoded elements from the stream into `out`.
// Returns the number of elements actually written (0 once exhausted).
inline int read_batch(StreamCtx& ctx, uint32_t* out, int max_count) {
    if (ctx.pos >= ctx.count || max_count <= 0) return 0;
    int avail = static_cast<int>(ctx.count - ctx.pos);
    int n = (avail < max_count) ? avail : max_count;
    std::memcpy(out, ctx.decoded.data() + ctx.pos, n * sizeof(uint32_t));
    ctx.pos += static_cast<uint32_t>(n);
    return n;
}

inline bool stream_has_more(const StreamCtx& ctx) {
    return ctx.pos < ctx.count;
}

// === .kpx candidate-set-driven decode (v7) ===
//
// Given a sorted candidate seq_id array, decode the .kpx posting and return
// the position vector for each candidate.  out is resized to candidates.size();
// out[i] holds the positions of candidates[i] (empty if the candidate is not
// present in the posting).  Returns false on corrupt input.
bool open_stream_kpx_for_candidates(
    const uint8_t* posting, size_t bytes,
    const uint32_t* candidates, size_t n_candidates,
    std::vector<std::vector<uint32_t>>& out);

// === posting blob inspection (no decode) ===

// Read the distinct_count u32 at the start of a .kix posting blob.  Returns
// 0 if the blob is shorter than the 8-byte header.
inline uint32_t posting_count(const uint8_t* posting, size_t bytes) {
    if (bytes < kPostingHeaderBytes) return 0;
    uint32_t cnt;
    std::memcpy(&cnt, posting, sizeof(uint32_t));
    return cnt;
}

// Name of the FastPFor ISA tier that the runtime dispatcher selected
// ("sse42" / "avx2" / "avx512bw" / "avx512vbmi2"; Phase 5f 4-tier ladder).
// First call resolves the tier; subsequent calls are cached.  Useful for
// diagnostic logging.
const char* active_tier_name();

} // namespace ikafssn::pfd
