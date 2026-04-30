#pragma once

#include <cstdint>
#include <cstring>
#include <utility>
#include <vector>
#include "index/pfd_codec.hpp"

namespace ikafssn {

// Streaming decoder for v7 .kpx postings (candidate-set-driven).
//
// .kpx v7 stores absolute positions partitioned per-(kmer, seq_id) above
// a build-time threshold, plus a self-describing short bucket carrying
// its own seq_id list and per-seq_id occurrence counts.  The decoder
// takes a sorted candidate seq_id array and returns per-candidate
// position vectors (out[i] holds positions for candidates[i], or an
// empty vector if candidates[i] does not appear in the posting).
//
// Construction takes a byte range that fully contains the k-mer's
// posting blob and a candidate seq_id array.  Decoding is deferred until
// the first call to positions_for() / for_each_candidate().
class PosDecoder {
public:
    PosDecoder() = default;
    PosDecoder(const uint8_t* data, const uint8_t* end,
               const uint32_t* candidates, std::size_t n_candidates)
        : data_(data),
          bytes_(end && end >= data ? static_cast<std::size_t>(end - data) : 0),
          candidates_(candidates),
          n_candidates_(n_candidates) {}

    // Materialise per-candidate position vectors (idempotent).
    void ensure_decoded() {
        if (decoded_) return;
        decoded_ = true;
        if (data_ && bytes_ > 0 && n_candidates_ > 0) {
            pfd::open_stream_kpx_for_candidates(data_, bytes_,
                                                candidates_, n_candidates_,
                                                per_cand_);
        } else {
            per_cand_.assign(n_candidates_, std::vector<uint32_t>{});
        }
    }

    // Position vector for a candidate index (0 <= i < n_candidates).
    const std::vector<uint32_t>& positions_for(std::size_t i) {
        ensure_decoded();
        return per_cand_[i];
    }

    // Iterate over all candidates that have at least one position.
    // cb is invoked as cb(seq_id, positions[]).
    template <class CB>
    void for_each_candidate(CB&& cb) {
        ensure_decoded();
        for (std::size_t i = 0; i < n_candidates_; i++) {
            if (!per_cand_[i].empty()) {
                cb(candidates_[i], per_cand_[i]);
            }
        }
    }

private:
    const uint8_t* data_ = nullptr;
    std::size_t bytes_ = 0;
    const uint32_t* candidates_ = nullptr;
    std::size_t n_candidates_ = 0;
    std::vector<std::vector<uint32_t>> per_cand_;
    bool decoded_ = false;
};

} // namespace ikafssn
