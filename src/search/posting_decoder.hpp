#pragma once

#include <cstdint>
#include "core/varint.hpp"

namespace ikafssn {

// Streaming decoder for delta-compressed position postings.
// Must be used in lockstep with SeqIdDecoder: call next(was_new_seq)
// where was_new_seq comes from the corresponding SeqIdDecoder.
class PosDecoder {
public:
    PosDecoder() = default;
    explicit PosDecoder(const uint8_t* data) : ptr_(data) {}

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

private:
    const uint8_t* ptr_ = nullptr;
    uint32_t prev_pos_ = 0;
};

} // namespace ikafssn
