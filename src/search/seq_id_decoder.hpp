#pragma once

#include <algorithm>
#include <cstdint>
#include "index/pfd_codec.hpp"

namespace ikafssn {

// Streaming decoder for v6 .kix postings.
//
// As of Phase 5e the on-disk format is custom (v5/v6: PForDelta on .kix,
// v6: per-(kmer, seq_id) partition + short bucket on .kpx).  The codec
// already returns absolute seq_ids in StreamCtx::decoded.  This decoder is
// therefore a thin "absolute id iterator" — it does not re-accumulate
// deltas.
//
// `was_new_seq()` reflects whether the last emitted id differs from the
// previous emitted id (a sequence boundary in the posting stream); it is
// preserved here for source compatibility with callers that still consume
// it, even though Phase 5c removed PosDecoder's dependence on it.
//
// Phase 5g-2: PosDecoder needs the decoded seq_id array up front so its
// .kpx merge can walk lock-step with the .kix stream.  ensure_decoded()
// is therefore exposed and decoded_data() / decoded_count() return the
// already-materialised buffer (zero copy).
class SeqIdDecoder {
public:
    static constexpr int kMaxBatch = 16;

    SeqIdDecoder() = default;
    explicit SeqIdDecoder(const uint8_t* data) : data_(data), bytes_(0) {}
    SeqIdDecoder(const uint8_t* data, const uint8_t* end)
        : data_(data),
          bytes_(end && end >= data ? static_cast<std::size_t>(end - data) : 0) {}

    // Force the underlying posting blob to be decoded into ctx_.decoded.
    // Idempotent; safe to call before any of the iterator methods or
    // before invoking decoded_data() / decoded_count().
    void ensure_decoded() {
        if (decoded_) return;
        decoded_ = true;
        if (data_ && bytes_ > 0) {
            pfd::open_stream_kix(data_, bytes_, ctx_);
        }
    }

    bool has_more() {
        ensure_decoded();
        return ctx_.pos < ctx_.count;
    }

    // Decode next absolute seq_id.
    uint32_t next() {
        ensure_decoded();
        uint32_t id = 0;
        if (ctx_.pos < ctx_.count) {
            id = ctx_.decoded[ctx_.pos++];
        }
        if (first_) {
            prev_id_ = id;
            first_ = false;
            was_new_seq_ = true;
        } else {
            was_new_seq_ = (id != prev_id_);
            prev_id_ = id;
        }
        return id;
    }

    bool was_new_seq() const { return was_new_seq_; }

    // For v4+, ptr() is no longer meaningful (data is pre-decoded).
    // Retained as a stub for source compatibility.
    const uint8_t* ptr() const { return data_; }

    // Phase 5g-2 accessors: zero-copy view of the decoded seq_id array.
    // Caller must invoke ensure_decoded() (or any iterator method) first.
    const uint32_t* decoded_data() const { return ctx_.decoded.data(); }
    std::size_t decoded_count() const { return ctx_.count; }

    // Decode up to max_count absolute seq_ids.  out_was_new[i] = 1 if the
    // emitted id differs from the previous emitted id (or this is the
    // first id overall).
    int next_batch(uint32_t* out_sids, uint8_t* out_was_new, int max_count) {
        ensure_decoded();
        if (max_count <= 0 || ctx_.pos >= ctx_.count) return 0;
        if (max_count > kMaxBatch) max_count = kMaxBatch;
        int avail = static_cast<int>(ctx_.count - ctx_.pos);
        int n = (avail < max_count) ? avail : max_count;
        for (int i = 0; i < n; i++) {
            uint32_t id = ctx_.decoded[ctx_.pos + i];
            if (first_) {
                prev_id_ = id;
                first_ = false;
                out_was_new[i] = 1;
            } else {
                out_was_new[i] = (id != prev_id_) ? 1 : 0;
                prev_id_ = id;
            }
            out_sids[i] = id;
        }
        ctx_.pos += n;
        if (n > 0) was_new_seq_ = (out_was_new[n - 1] != 0);
        return n;
    }

private:
    const uint8_t* data_ = nullptr;
    std::size_t bytes_ = 0;
    pfd::StreamCtx ctx_;
    bool decoded_ = false;
    uint32_t prev_id_ = 0;
    bool first_ = true;
    bool was_new_seq_ = false;
};

} // namespace ikafssn
