#pragma once

// Phase 5g-1 — split codec for .kix v5 and .kpx v5.
//
//   .kix v5: FastPFor CompositeCodec<SIMDFastPFor<4>, VariableByte>
//            (PForDelta + VByte tail).  Restored from Phase 5b after the
//            Phase 5e custom block codec was found to mishandle the
//            outlier-in-block case that PForDelta's exception stream is
//            specifically designed for.  Posting layout on disk:
//              [u32 count]                  — number of seq_id deltas
//              [u32 payload_words]          — payload size in u32 words
//              [u32 payload[payload_words]] — codec output, byte-unaligned
//                                             on the wire (use memcpy)
//
//   .kpx v5: custom FOR-within-block codec (Phase 5e), unchanged in v5.
//            Each 128-element block subtracts its min before bitpacking,
//            making the per-block bit-width depend on within-sequence
//            spread rather than absolute position magnitude.  Posting
//            layout on disk:
//              [u32 count]
//              repeated count/128 times:
//                [u8 b][u32 min][128*b/8 bytes bitpacked (value-min)]
//              [u8 tail_count][u32 tail_min][varint stream of tail-min]
//
// Both reader entry points materialise the entire decoded posting into
// StreamCtx::decoded at open_stream_*() time; refill() is a no-op kept
// for API symmetry with potential future per-block streaming codecs.

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

// Name of the FastPFor ISA tier that the runtime dispatcher selected
// ("sse42" / "avx2" / "avx512bw" / "avx512vbmi2"; Phase 5f 4-tier ladder).
// First call resolves the tier; subsequent calls are cached.  Useful for
// diagnostic logging.
const char* active_tier_name();

} // namespace ikafssn::pfd
