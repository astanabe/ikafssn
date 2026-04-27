#pragma once

#include <algorithm>
#include <cstdint>
#include "core/varint.hpp"
#include "index/varint_simd.hpp"

namespace ikafssn {

// Streaming decoder for delta-compressed ID postings.
// Decodes one seq_id at a time (next()) or up to N at a time (next_batch()).
class SeqIdDecoder {
public:
    static constexpr int kMaxBatch = kVarintMaxBatch;

    SeqIdDecoder() = default;
    explicit SeqIdDecoder(const uint8_t* data) : ptr_(data), end_(nullptr) {}
    SeqIdDecoder(const uint8_t* data, const uint8_t* end) : ptr_(data), end_(end) {}

    bool has_more() const { return end_ && ptr_ < end_; }

    // Decode next seq_id. Returns the absolute seq_id.
    uint32_t next() {
        uint32_t delta;
        ptr_ += varint_decode(ptr_, delta);
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

    // Did the last next()/next_batch() decode a new (different) seq_id?
    // For next_batch, reflects the last element in the returned batch.
    // Used by PosDecoder to detect sequence boundaries for delta reset.
    bool was_new_seq() const { return was_new_seq_; }

    // Current position in the byte stream.
    const uint8_t* ptr() const { return ptr_; }

    // Decode up to max_count seq_ids in a single batched call.
    //   out_sids[i]    = absolute seq_id of the i-th decoded element
    //   out_was_new[i] = 1 if delta != 0 (or first element overall), else 0
    // Returns the number of elements actually decoded (0..max_count). Updates
    // prev_id_ / first_ / was_new_seq_ as a side-effect (last applied value).
    int next_batch(uint32_t* out_sids, uint8_t* out_was_new, int max_count) {
        if (!end_ || ptr_ >= end_ || max_count <= 0) return 0;
        if (max_count > kMaxBatch) max_count = kMaxBatch;
        uint32_t deltas[kMaxBatch];
        std::size_t consumed = 0;
        int n = varint_decode_batch(ptr_, end_, deltas, &consumed, max_count);
        ptr_ += consumed;

        for (int i = 0; i < n; i++) {
            uint32_t delta = deltas[i];
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
        if (n > 0) was_new_seq_ = (out_was_new[n - 1] != 0);
        return n;
    }

private:
    const uint8_t* ptr_ = nullptr;
    const uint8_t* end_ = nullptr;
    uint32_t prev_id_ = 0;
    bool first_ = true;
    bool was_new_seq_ = false;
};

} // namespace ikafssn
