#pragma once

#include <cstdint>
#include "core/varint.hpp"

namespace ikafssn {

// Streaming decoder for delta-compressed ID postings.
// Decodes one seq_id at a time.
class SeqIdDecoder {
public:
    SeqIdDecoder() = default;
    explicit SeqIdDecoder(const uint8_t* data) : ptr_(data) {}

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

    // Did the last next() decode a new (different) seq_id?
    // Used by PosDecoder to detect sequence boundaries for delta reset.
    bool was_new_seq() const { return was_new_seq_; }

    // Current position in the byte stream.
    const uint8_t* ptr() const { return ptr_; }

private:
    const uint8_t* ptr_ = nullptr;
    uint32_t prev_id_ = 0;
    bool first_ = true;
    bool was_new_seq_ = false;
};

} // namespace ikafssn
