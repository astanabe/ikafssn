#pragma once

// Phase 5b — SIMD-FastPFOR* + Simple-8b/VByte posting codec (.kix v4 / .kpx v4).
//
// Posting layout (per k-mer, on disk):
//
//   [u32 count]              — number of postings encoded
//   [u32 payload_words]      — number of uint32_t words in the encoded payload
//   [uint32 payload[N]]      — encoded payload (FastPFor CompositeCodec output)
//
// For .kix the input array is the **delta-encoded** seq_id stream (first
// element absolute, then differences).  For .kpx the input array is the
// **absolute** position stream — FastPFor's per-block bit-width adapts to
// the position range so dense per-sequence clusters compress efficiently
// without the explicit within-seq delta reset that v3 needed.
//
// The on-disk byte stream is unaligned wrt 4-byte boundaries; readers use
// std::memcpy to materialise the count + payload-word-count + payload
// pointer.  Decoded data is materialised lazily into the StreamCtx buffer
// when open_stream_*() is called.

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

// Per-posting byte-stream header (count + payload_words). Posting begins
// with these two u32 values written little-endian; payload follows.
inline constexpr size_t kPostingHeaderBytes = 8;

// === posting-level encode wrappers ===

// Encode the seq_id-delta stream for a .kix posting. Writes count +
// payload_words + payload into `out` (appended).  Returns the number of
// bytes written.
size_t encode_posting_kix(const uint32_t* delta_array, uint32_t count,
                          std::vector<uint8_t>& out);

// Encode the absolute-position stream for a .kpx posting. Writes count +
// payload_words + payload into `out` (appended).  Returns bytes written.
size_t encode_posting_kpx(const uint32_t* abs_pos_array, uint32_t count,
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

// Initialise a StreamCtx for a .kpx posting.  Symmetric to open_stream_kix.
bool open_stream_kpx(const uint8_t* posting, size_t bytes, StreamCtx& ctx);

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

} // namespace ikafssn::pfd
