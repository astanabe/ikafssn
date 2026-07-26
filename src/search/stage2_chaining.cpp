#include "search/stage2_chaining.hpp"
#include "search/diagonal_filter.hpp"

#include <algorithm>
#include <cstdlib>

namespace ikafssn {

namespace {

// Per-thread working buffers for chain_hits().  They grow monotonically to the
// largest subject seen so far, so the steady state performs no allocation.
struct ChainScratch {
    std::vector<Hit>      work;           // sorted, deduped, filtered hits
    std::vector<uint32_t> q;              // DP input, split out of `work`
    std::vector<uint32_t> s;
    std::vector<int32_t>  diag;           // s - q, per hit
    std::vector<uint32_t> dp;
    std::vector<int32_t>  prev;
    std::vector<size_t>   chain_indices;  // traceback, strictly decreasing
    std::vector<int32_t>  diag_counts;    // diagonal_filter counting buffer
};

inline ChainScratch& tls_chain_scratch() {
    thread_local ChainScratch scratch;
    return scratch;
}

// Make `v` hold at least n elements and hand back its data pointer.  Callers
// overwrite every element they read, so the array never shrinks and the grown
// part needs no particular value.
template <typename T>
inline T* scratch_ptr(std::vector<T>& v, size_t n) {
    if (v.size() < n) v.resize(n);
    return v.data();
}

// Function objects rather than free functions: std::sort / std::unique would
// otherwise be instantiated on a function pointer and call through it.
struct HitLess {
    bool operator()(const Hit& a, const Hit& b) const {
        return a.q_pos < b.q_pos || (a.q_pos == b.q_pos && a.s_pos < b.s_pos);
    }
};

struct HitSamePos {
    bool operator()(const Hit& a, const Hit& b) const {
        return a.q_pos == b.q_pos && a.s_pos == b.s_pos;
    }
};

}  // namespace

void chain_hits(const std::vector<Hit>& raw_hits,
                SeqId seq_id,
                int span,
                bool is_reverse,
                const Stage2Config& config,
                std::vector<ChainResult>& out) {
    out.clear();
    if (raw_hits.empty()) return;

    ChainScratch& scratch = tls_chain_scratch();
    std::vector<Hit>& work = scratch.work;
    std::vector<size_t>& chain_indices = scratch.chain_indices;

    // Deduplicate (q_pos, s_pos) pairs from degenerate base expansion.
    work.assign(raw_hits.begin(), raw_hits.end());
    std::sort(work.begin(), work.end(), HitLess{});
    work.erase(std::unique(work.begin(), work.end(), HitSamePos{}), work.end());

    diagonal_filter(work, config.min_nhit_diag, scratch.diag_counts);
    if (work.empty()) return;

    // Extraction limit and tie behaviour.  N == 0 means unlimited.  In a
    // tie-inclusive mode the chainscore of the N-th chain (s_n) bounds how
    // far past N we keep extracting: chainscore is non-increasing across
    // iterations, so once a chain below s_n appears we stop.
    const uint32_t n_limit = config.max_nhit_per_subject;
    const bool unlimited = (n_limit == 0);
    const bool tie_inclusive = (config.max_nhit_per_subject_mode == 3 ||
                                config.max_nhit_per_subject_mode == 4);
    uint32_t s_n = 0;
    bool s_n_set = false;

    while (!work.empty()) {
        // Strict take-N stop: N chains already produced.
        if (!unlimited && !tie_inclusive && out.size() >= n_limit) break;

        const size_t n = work.size();

        // Split the hits into parallel arrays and precompute each hit's
        // diagonal.  The gap test between hits i and j is
        //   |(s_i - s_j) - (q_i - q_j)| = |diag_i - diag_j|,
        // so the DP inner loop reduces to one subtraction on values it can
        // read straight out of `diag`.  The buffers live in thread-local
        // storage, so index them through local pointers to keep the loop free
        // of base-pointer reloads.
        const Hit* const wp = work.data();
        uint32_t* const qp = scratch_ptr(scratch.q, n);
        uint32_t* const sp = scratch_ptr(scratch.s, n);
        int32_t* const dg = scratch_ptr(scratch.diag, n);
        uint32_t* const dpp = scratch_ptr(scratch.dp, n);
        int32_t* const prevp = scratch_ptr(scratch.prev, n);
        for (size_t i = 0; i < n; i++) {
            qp[i] = wp[i].q_pos;
            sp[i] = wp[i].s_pos;
            dg[i] = static_cast<int32_t>(wp[i].s_pos) -
                    static_cast<int32_t>(wp[i].q_pos);
        }

        const int64_t max_gap = static_cast<int64_t>(config.max_gap);
        const size_t lookback = config.chain_max_lookback;

        dpp[0] = 1;
        prevp[0] = -1;
        for (size_t i = 1; i < n; i++) {
            const size_t j_start = (lookback > 0 && i > lookback) ? (i - lookback) : 0;
            const uint32_t q_i = qp[i];
            const uint32_t s_i = sp[i];
            const int64_t diag_i = dg[i];
            // Predecessors are scanned in ascending j and only a strictly
            // longer chain replaces the incumbent, so ties keep the smallest
            // j.  Reversing or vectorising this scan would change which
            // predecessor wins a tie, and with it the reported coordinates.
            uint32_t dp_i = 1;
            int32_t prev_i = -1;
            for (size_t j = j_start; j < i; j++) {
                if (qp[j] >= q_i) continue;
                if (sp[j] >= s_i) continue;
                if (dpp[j] < dp_i) continue;
                if (std::llabs(diag_i - static_cast<int64_t>(dg[j])) > max_gap) continue;
                dp_i = dpp[j] + 1;
                prev_i = static_cast<int32_t>(j);
            }
            dpp[i] = dp_i;
            prevp[i] = prev_i;
        }

        // Find best chain endpoint
        size_t best_idx = 0;
        for (size_t i = 1; i < n; i++) {
            if (dpp[i] > dpp[best_idx]) {
                best_idx = i;
            }
        }

        uint32_t best_score = dpp[best_idx];
        if (best_score < config.min_score) break;

        // Tie-inclusive stop: past the N-th chain, only chains tying the
        // N-th chainscore survive.
        if (!unlimited && tie_inclusive && s_n_set && best_score < s_n) break;

        // Traceback to collect chain hit indices
        chain_indices.clear();
        size_t idx = best_idx;
        while (idx != SIZE_MAX) {
            chain_indices.push_back(idx);
            idx = (prevp[idx] >= 0) ? static_cast<size_t>(prevp[idx]) : SIZE_MAX;
        }

        // chain_indices is in reverse order; first element is end, last is start
        size_t chain_start_idx = chain_indices.back();
        size_t chain_end_idx = chain_indices.front();

        ChainResult cr{};
        cr.seq_id = seq_id;
        cr.chainscore = best_score;
        cr.is_reverse = is_reverse;
        cr.q_start = work[chain_start_idx].q_pos;
        cr.q_end = work[chain_end_idx].q_pos + static_cast<uint32_t>(span);
        cr.s_start = work[chain_start_idx].s_pos;
        cr.s_end = work[chain_end_idx].s_pos + static_cast<uint32_t>(span);

        out.push_back(cr);

        // Record the N-th chain's chainscore the moment N is reached
        // (tie-inclusive mode only).
        if (!unlimited && tie_inclusive && !s_n_set &&
            out.size() == n_limit) {
            s_n = best_score;
            s_n_set = true;
        }

        // Early return for a strict single-chain take (no removal overhead).
        if (!unlimited && !tie_inclusive && n_limit == 1) break;

        // Drop the chain's hits so the next iteration sees the remainder.
        // prev[i] < i holds for every link, so chain_indices is strictly
        // decreasing and walking it backwards yields the removal targets in
        // ascending order — enough to compact `work` in a single pass.
        size_t w = 0;
        size_t ci = chain_indices.size();
        for (size_t r = 0; r < n; r++) {
            if (ci > 0 && chain_indices[ci - 1] == r) { ci--; continue; }
            work[w++] = work[r];
        }
        work.resize(w);
    }
}

} // namespace ikafssn
