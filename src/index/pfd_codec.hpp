#pragma once

// Phase 5g-2 — split codec for .kix v6 and .kpx v6.
//
//   .kix v6: FastPFor CompositeCodec<SIMDFastPFor<4>, VariableByte>
//            (PForDelta + VByte tail).  Unchanged on the wire from the
//            Phase 5g-1 v5 layout; bumped to v6 because the rest of the
//            index family (Phase 5g-2 .kpx layout change) bumps too.
//            Posting layout on disk:
//              [u32 count]                  — number of seq_id deltas
//              [u32 payload_words]          — payload size in u32 words
//              [u32 payload[payload_words]] — codec output, byte-unaligned
//                                             on the wire (use memcpy)
//
//   .kpx v6: per-(kmer, seq_id) partitioned FOR-within-block layout.
//            Each (k-mer, seq_id) cluster whose occurrence count exceeds
//            the build-time `freq_threshold_part` is split out into its
//            own partition group; remaining occurrences are merged into a
//            single short-bucket FOR-block stream that absorbs the long
//            tail of low-multiplicity (k-mer, seq_id) cells.  This breaks
//            the coupling between absolute position magnitude and the
//            per-block bit-width that hurt the Phase 5e v5 layout on
//            chromosome-class subjects.  Posting layout on disk:
//              [u32 total_count]
//              [u32 partition_count]
//              repeated partition_count times (sorted by seq_id):
//                [u32 seq_id]
//                [u32 occurrence_count]
//                [FOR-block stream over occurrence_count positions]
//              [u32 short_bucket_count]
//              [FOR-block stream over short_bucket_count positions]
//            Each FOR-block stream is the same encode_block_for + tail
//            FOR layout used since Phase 5e:
//                repeated count/128 times:
//                  [u8 b][u32 min][128*b/8 bytes bitpacked (value-min)]
//                [u8 tail_count][u32 tail_min][varint stream of tail-min]
//            (tail_min/varint are emitted only when tail_count > 0.)
//
// .kix open_stream materialises the entire decoded posting (absolute
// seq_ids) into StreamCtx::decoded.  .kpx open_stream takes the caller-
// provided sid_stream (typically the .kix decoded array for the same
// k-mer) and emits positions in lock-step with that stream — the partition
// groups are scanned monotonically with a per-group cursor and the short
// bucket fills the gaps.

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

// Per-posting byte-stream header for .kix (count + payload_words). The
// .kpx layout uses its own variable-length header (see the file-level
// comment above) so this constant only describes the .kix posting blob.
inline constexpr size_t kPostingHeaderBytes = 8;

// === posting-level encode wrappers ===

// Encode the seq_id-delta stream for a .kix posting. Writes count +
// payload_words + payload into `out` (appended).  Returns the number of
// bytes written.
size_t encode_posting_kix(const uint32_t* delta_array, uint32_t count,
                          std::vector<uint8_t>& out);

// Encode the absolute-position stream for a .kpx posting (v6).
//   sid_array:               length = count, sorted (may have duplicates)
//   abs_pos_array:           length = count, parallel to sid_array
//   freq_threshold_part:   per-(kmer, seq_id) count > threshold => the
//                            cluster becomes its own partition group;
//                            otherwise its positions go into the short
//                            bucket
// The sid_array merely partitions the positions on the encode side; only
// the position stream is stored (the merge against sid_array happens at
// decode time, driven by the caller-provided .kix seq_id stream).
size_t encode_posting_kpx(const uint32_t* sid_array,
                          const uint32_t* abs_pos_array,
                          uint32_t count,
                          uint32_t freq_threshold_part,
                          std::vector<uint8_t>& out);

// === streaming decode context ===
//
// open_stream_kix/_kpx materialises the entire decoded posting into
// `decoded`.  refill() is a no-op for the wrapper backend (the data is
// already fully decoded) — it's retained for API symmetry with future
// per-block streaming codecs.
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

// Initialise a StreamCtx for a .kpx posting (v6).  The decoded stream is
// merged in lock-step with the caller-supplied `sid_stream` (typically the
// .kix-decoded seq_id array for the same k-mer): on return,
// ctx.decoded[i] is the position corresponding to sid_stream[i].  Returns
// false on header/size mismatch or if the partition-group seq_ids do not
// align with `sid_stream` (corrupt index).
bool open_stream_kpx(const uint8_t* posting, size_t bytes,
                     const uint32_t* sid_stream, size_t n_sids,
                     StreamCtx& ctx);

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

// === posting blob inspection (no decode) ===

// Read the count u32 at the start of a posting blob.  Returns 0 if the
// blob is shorter than the 8-byte header.
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
