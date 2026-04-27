#pragma once

#include <algorithm>
#include <cstdint>
#include "core/varint.hpp"
#include "index/varint_simd.hpp"

namespace ikafssn {

// Streaming decoder for delta-compressed position postings.
// Must be used in lockstep with SeqIdDecoder: call next(was_new_seq) where
// was_new_seq comes from the corresponding SeqIdDecoder, or use next_batch()
// with the matching was_new_seq[] array from SeqIdDecoder::next_batch().
class PosDecoder {
public:
    static constexpr int kMaxBatch = kVarintMaxBatch;

    PosDecoder() = default;
    explicit PosDecoder(const uint8_t* data) : ptr_(data), end_(nullptr) {}
    PosDecoder(const uint8_t* data, const uint8_t* end) : ptr_(data), end_(end) {}

    bool has_more() const { return end_ && ptr_ < end_; }

    // Decode next position. was_new_seq indicates sequence boundary (delta reset).
    uint32_t next(bool was_new_seq) {
        uint32_t val;
        ptr_ += varint_decode(ptr_, val);
        if (was_new_seq) {
            prev_pos_ = val; // raw value, not delta
        } else {
            prev_pos_ += val; // delta from previous pos
        }
        return prev_pos_;
    }

    const uint8_t* ptr() const { return ptr_; }

    // Decode up to max_count positions. was_new_seq[i] (typically passed from
    // the corresponding SeqIdDecoder::next_batch() output) controls per-element
    // delta reset. Returns number decoded.
    //
    // When constructed with an end pointer, uses the SIMD-friendly
    // varint_decode_batch path. When end_ is null (legacy construction),
    // falls back to a pure scalar varint_decode loop bounded only by
    // max_count, matching the original scalar next() exactly. The lockstep
    // SeqIdDecoder is responsible for keeping max_count within the producer's
    // length so we never read past the legitimate stream.
    int next_batch(uint32_t* out_pos, const uint8_t* was_new_seq, int max_count) {
        if (max_count <= 0) return 0;
        if (max_count > kMaxBatch) max_count = kMaxBatch;

        if (end_) {
            uint32_t raw[kMaxBatch];
            std::size_t consumed = 0;
            int n = varint_decode_batch(ptr_, end_, raw, &consumed, max_count);
            ptr_ += consumed;
            for (int i = 0; i < n; i++) {
                if (was_new_seq[i]) {
                    prev_pos_ = raw[i];
                } else {
                    prev_pos_ += raw[i];
                }
                out_pos[i] = prev_pos_;
            }
            return n;
        }

        // Unbounded: scalar decode, bit-exact with original next() x max_count.
        for (int i = 0; i < max_count; i++) {
            uint32_t v;
            ptr_ += varint_decode(ptr_, v);
            if (was_new_seq[i]) {
                prev_pos_ = v;
            } else {
                prev_pos_ += v;
            }
            out_pos[i] = prev_pos_;
        }
        return max_count;
    }

private:
    const uint8_t* ptr_ = nullptr;
    const uint8_t* end_ = nullptr;
    uint32_t prev_pos_ = 0;
};

} // namespace ikafssn
