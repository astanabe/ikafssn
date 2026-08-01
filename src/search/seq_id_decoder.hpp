#pragma once

#include <algorithm>
#include <cstdint>
#include "index/pfd_codec.hpp"

namespace ikafssn {

// Decoder for one .kix posting list.
//
// A .kix posting list holds **distinct** seq_ids only (intra-sequence k-mer
// duplicates are removed by a SIMD dedup kernel at build time).  This class
// owns the decode buffer and exposes the decoded absolute seq_ids as a
// zero-copy span.
class SeqIdDecoder {
public:
    SeqIdDecoder() = default;
    SeqIdDecoder(const uint8_t* data, const uint8_t* end)
        : data_(data),
          bytes_(end && end >= data ? static_cast<std::size_t>(end - data) : 0) {}

    // Point at a new posting list without reconstructing the object.  The
    // ctx_.decoded buffer is left intact for the next open_stream_kix to
    // decode into.
    void reset(const uint8_t* data, const uint8_t* end) {
        data_ = data;
        bytes_ = (end && end >= data)
                     ? static_cast<std::size_t>(end - data)
                     : 0;
        decoded_   = false;
        ctx_.count = 0;
    }

    void ensure_decoded() {
        if (decoded_) return;
        decoded_ = true;
        if (data_ && bytes_ > 0) {
            pfd::open_stream_kix(data_, bytes_, ctx_);
        }
    }

    // Zero-copy view of the decoded distinct seq_id array.  Valid until the
    // next reset() + ensure_decoded().
    const uint32_t* decoded_data() const { return ctx_.decoded.data(); }
    std::size_t decoded_count() const { return ctx_.count; }

private:
    const uint8_t* data_ = nullptr;
    std::size_t bytes_ = 0;
    pfd::StreamCtx ctx_;
    bool decoded_ = false;
};

} // namespace ikafssn
