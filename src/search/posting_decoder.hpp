#pragma once

#include <cstdint>
#include <cstring>
#include <utility>
#include <vector>
#include "index/pfd_codec.hpp"

namespace ikafssn {

// Streaming decoder for .kpx posting lists (candidate-set-driven).
//
// .kpx stores absolute positions partitioned per-(kmer, seq_id) above
// a build-time threshold and splits the short bucket into occ=1 and
// occ>=2 sub-buckets.  No seq_ids are stored in the .kpx posting list
// itself — the .kix decoded distinct seq_id array supplies the
// resolution between position rank and seq_id.  The decoder therefore
// needs:
//
//   - the .kix decoded distinct seq_id array (kix_decoded[0..kix_count),
//     strictly increasing — provided by SeqIdDecoder::decoded_data() /
//     decoded_count())
//   - a sorted candidate seq_id array (subset of kix_decoded)
//   - a per-thread reusable PosDecodeScratch
//
// All three pointers must outlive the PosDecoder.  Decoding is deferred
// until the first call to positions_for() / for_each_candidate().
class PosDecoder {
public:
    PosDecoder() = default;
    PosDecoder(const uint8_t* data, const uint8_t* end,
               const uint32_t* kix_decoded, std::size_t kix_count,
               const uint32_t* candidates, std::size_t n_candidates,
               pfd::PosDecodeScratch* scratch)
        : data_(data),
          bytes_(end && end >= data ? static_cast<std::size_t>(end - data) : 0),
          kix_decoded_(kix_decoded),
          kix_count_(kix_count),
          candidates_(candidates),
          n_candidates_(n_candidates),
          scratch_(scratch) {}

    // Materialise per-candidate position vectors (idempotent).
    void ensure_decoded() {
        if (decoded_) return;
        decoded_ = true;
        // kix_count_ == 0 means the k-mer has no .kix posting list, so its
        // .kpx posting list cannot be addressed at all (its pos_offset may
        // be 0 — aliasing the first k-mer's posting list — because the
        // builder does not write placeholders for empty k-mers).  Bail
        // out before touching the .kpx data.
        if (data_ && bytes_ > 0 && n_candidates_ > 0 && kix_count_ > 0 && scratch_) {
            pfd::open_stream_kpx_for_candidates(data_, bytes_,
                                                kix_decoded_, kix_count_,
                                                candidates_, n_candidates_,
                                                *scratch_,
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
    const uint32_t* kix_decoded_ = nullptr;
    std::size_t kix_count_ = 0;
    const uint32_t* candidates_ = nullptr;
    std::size_t n_candidates_ = 0;
    pfd::PosDecodeScratch* scratch_ = nullptr;
    std::vector<std::vector<uint32_t>> per_cand_;
    bool decoded_ = false;
};

} // namespace ikafssn
