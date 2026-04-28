#pragma once

#include <algorithm>
#include <cstdint>
#include <cstring>
#include "index/pfd_codec.hpp"

namespace ikafssn {

// Streaming decoder for v6 .kpx postings.
//
// .kpx v6 stores absolute positions partitioned per-(kmer, seq_id) above
// a build-time threshold, plus a shared short bucket for the rest (see
// src/index/pfd_codec.hpp).  The decoder requires the caller to provide
// the .kix-decoded seq_id array for the same k-mer so the merge can walk
// lock-step against it; ctx_.decoded[i] then corresponds to sids[i].
//
// Construction takes a byte range that fully contains the k-mer's
// posting blob (the leading [u32 total_count][u32 partition_count]
// header — and the per-section trailers — describes the actual blob
// length).
class PosDecoder {
public:
    static constexpr int kMaxBatch = 16;

    PosDecoder() = default;
    PosDecoder(const uint8_t* data, const uint8_t* end,
               const uint32_t* sids, std::size_t n_sids)
        : data_(data),
          bytes_(end && end >= data ? static_cast<std::size_t>(end - data) : 0),
          sids_(sids),
          n_sids_(n_sids) {}

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
            pfd::open_stream_kpx(data_, bytes_, sids_, n_sids_, ctx_);
        }
    }

    const uint8_t* data_ = nullptr;
    std::size_t bytes_ = 0;
    const uint32_t* sids_ = nullptr;
    std::size_t n_sids_ = 0;
    pfd::StreamCtx ctx_;
    bool decoded_ = false;
};

} // namespace ikafssn
