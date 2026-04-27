#pragma once

#include <algorithm>
#include <cstdint>
#include <cstring>
#include "index/pfd_codec.hpp"

namespace ikafssn {

// Streaming decoder for v4 .kpx postings.
//
// .kpx v4 stores absolute positions (FastPFor-encoded), so unlike v3 the
// decoder does not need any sequence-boundary signal from the
// SeqIdDecoder.  Construction takes a byte range that fully contains the
// k-mer's posting blob (the leading [u32 count][u32 payload_words] header
// describes the actual blob length).
class PosDecoder {
public:
    static constexpr int kMaxBatch = 16;

    PosDecoder() = default;
    PosDecoder(const uint8_t* data, const uint8_t* end)
        : data_(data),
          bytes_(end && end >= data ? static_cast<std::size_t>(end - data) : 0) {}

    bool has_more() {
        ensure_decoded();
        return ctx_.pos < ctx_.count;
    }

    uint32_t next() {
        ensure_decoded();
        if (ctx_.pos < ctx_.count) {
            return ctx_.decoded[ctx_.pos++];
        }
        return 0;
    }

    int next_batch(uint32_t* out_pos, int max_count) {
        ensure_decoded();
        if (max_count <= 0 || ctx_.pos >= ctx_.count) return 0;
        if (max_count > kMaxBatch) max_count = kMaxBatch;
        int avail = static_cast<int>(ctx_.count - ctx_.pos);
        int n = (avail < max_count) ? avail : max_count;
        std::memcpy(out_pos, ctx_.decoded.data() + ctx_.pos, n * sizeof(uint32_t));
        ctx_.pos += n;
        return n;
    }

private:
    void ensure_decoded() {
        if (decoded_) return;
        decoded_ = true;
        if (data_ && bytes_ > 0) {
            pfd::open_stream_kpx(data_, bytes_, ctx_);
        }
    }

    const uint8_t* data_ = nullptr;
    std::size_t bytes_ = 0;
    pfd::StreamCtx ctx_;
    bool decoded_ = false;
};

} // namespace ikafssn
