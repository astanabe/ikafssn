#pragma once

#include <algorithm>
#include <cstdint>
#include "index/pfd_codec.hpp"

namespace ikafssn {

// Streaming decoder for v4 .kix postings (FastPFor-encoded delta stream).
// API is preserved from the v3 implementation: next(), next_batch(), and
// was_new_seq() are all available, but the underlying storage is now a
// pre-decoded delta buffer (lazily materialised on first read) rather than
// a varint byte stream.  Construction takes the same `(data, end)` byte
// range so callers do not need to change.
class SeqIdDecoder {
public:
    static constexpr int kMaxBatch = 16;

    SeqIdDecoder() = default;
    explicit SeqIdDecoder(const uint8_t* data) : data_(data), bytes_(0) {}
    SeqIdDecoder(const uint8_t* data, const uint8_t* end)
        : data_(data),
          bytes_(end && end >= data ? static_cast<std::size_t>(end - data) : 0) {}

    bool has_more() {
        ensure_decoded();
        return ctx_.pos < ctx_.count;
    }

    // Decode next seq_id.
    uint32_t next() {
        ensure_decoded();
        uint32_t delta = 0;
        if (ctx_.pos < ctx_.count) {
            delta = ctx_.decoded[ctx_.pos++];
        }
        if (first_) {
            prev_id_ = delta;
            first_ = false;
            was_new_seq_ = true;
        } else {
            was_new_seq_ = (delta != 0);
            prev_id_ += delta;
        }
        return prev_id_;
    }

    bool was_new_seq() const { return was_new_seq_; }

    // For v4, ptr() is no longer meaningful (data is pre-decoded). Retained
    // as a stub for source compatibility — callers that diff'd against ptr()
    // are flagged at compile time if they rely on its identity.
    const uint8_t* ptr() const { return data_; }

    // Decode up to max_count seq_ids.  out_was_new[i] = 1 if delta != 0
    // (sequence boundary) or the very first element overall.
    int next_batch(uint32_t* out_sids, uint8_t* out_was_new, int max_count) {
        ensure_decoded();
        if (max_count <= 0 || ctx_.pos >= ctx_.count) return 0;
        if (max_count > kMaxBatch) max_count = kMaxBatch;
        int avail = static_cast<int>(ctx_.count - ctx_.pos);
        int n = (avail < max_count) ? avail : max_count;
        for (int i = 0; i < n; i++) {
            uint32_t delta = ctx_.decoded[ctx_.pos + i];
            if (first_) {
                prev_id_ = delta;
                first_ = false;
                out_was_new[i] = 1;
            } else {
                out_was_new[i] = (delta != 0) ? 1 : 0;
                prev_id_ += delta;
            }
            out_sids[i] = prev_id_;
        }
        ctx_.pos += n;
        if (n > 0) was_new_seq_ = (out_was_new[n - 1] != 0);
        return n;
    }

private:
    void ensure_decoded() {
        if (decoded_) return;
        decoded_ = true;
        if (data_ && bytes_ > 0) {
            pfd::open_stream_kix(data_, bytes_, ctx_);
        }
    }

    const uint8_t* data_ = nullptr;
    std::size_t bytes_ = 0;
    pfd::StreamCtx ctx_;
    bool decoded_ = false;
    uint32_t prev_id_ = 0;
    bool first_ = true;
    bool was_new_seq_ = false;
};

} // namespace ikafssn
