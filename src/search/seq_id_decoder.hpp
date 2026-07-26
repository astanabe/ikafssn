#pragma once

#include <algorithm>
#include <cstdint>
#include "index/pfd_codec.hpp"

namespace ikafssn {

// Decoder for one .kix posting list.
//
// The .kix posting list stream contains **distinct** seq_ids only
// (intra-sequence k-mer duplicates are removed by a SIMD dedup kernel at
// build time).  The codec decodes the whole posting list into
// StreamCtx::decoded as absolute seq_ids; this class owns that buffer and
// exposes it as a zero-copy span so callers can consume the postings
// without an intermediate copy.
class SeqIdDecoder {
public:
    SeqIdDecoder() = default;
    SeqIdDecoder(const uint8_t* data, const uint8_t* end)
        : data_(data),
          bytes_(end && end >= data ? static_cast<std::size_t>(end - data) : 0) {}

    // Reset to point at a new posting list without reconstructing the
    // object. Reuses the StreamCtx::decoded heap buffer so the next
    // ensure_decoded() does not have to re-grow it from zero.
    void reset(const uint8_t* data, const uint8_t* end) {
        data_ = data;
        bytes_ = (end && end >= data)
                     ? static_cast<std::size_t>(end - data)
                     : 0;
        decoded_     = false;
        ctx_.count   = 0;
        ctx_.pos     = 0;
        // Keep ctx_.decoded capacity intact so the upcoming
        // open_stream_kix can resize() in-place.
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
