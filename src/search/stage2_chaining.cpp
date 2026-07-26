#include "search/stage2_chaining.hpp"
#include "search/diagonal_filter.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>

namespace ikafssn {

namespace {

// Per-thread working buffers for chain_hits().  They grow monotonically to the
// largest subject seen so far, so the steady state performs no allocation.
struct ChainScratch {
    std::vector<Hit>      work;           // sorted, deduped, filtered hits
    std::vector<uint32_t> dp;
    std::vector<int32_t>  prev;
    std::vector<size_t>   chain_indices;  // traceback, strictly decreasing
    std::vector<int32_t>  diags;          // diagonal_filter counting buffer
};

inline ChainScratch& tls_chain_scratch() {
    thread_local ChainScratch scratch;
    return scratch;
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
    std::vector<uint32_t>& dp = scratch.dp;
    std::vector<int32_t>& prev = scratch.prev;
    std::vector<size_t>& chain_indices = scratch.chain_indices;

    // Deduplicate (q_pos, s_pos) pairs from degenerate base expansion.
    work.assign(raw_hits.begin(), raw_hits.end());
    std::sort(work.begin(), work.end(), HitLess{});
    work.erase(std::unique(work.begin(), work.end(), HitSamePos{}), work.end());

    diagonal_filter(work, config.min_nhit_diag, scratch.diags);
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

        dp.assign(n, 1);
        prev.assign(n, -1);

        // The buffers live in thread-local storage, so index them through
        // local pointers to keep the DP loop free of base-pointer reloads.
        const Hit* const wp = work.data();
        uint32_t* const dpp = dp.data();
        int32_t* const prevp = prev.data();
        const int64_t max_gap = static_cast<int64_t>(config.max_gap);
        const size_t lookback = config.chain_max_lookback;

        for (size_t i = 1; i < n; i++) {
            size_t j_start = (lookback > 0 && i > lookback) ? (i - lookback) : 0;
            for (size_t j = j_start; j < i; j++) {
                if (wp[j].q_pos >= wp[i].q_pos) continue;
                if (wp[j].s_pos >= wp[i].s_pos) continue;

                int64_t gap_q = static_cast<int64_t>(wp[i].q_pos) - static_cast<int64_t>(wp[j].q_pos);
                int64_t gap_s = static_cast<int64_t>(wp[i].s_pos) - static_cast<int64_t>(wp[j].s_pos);
                int64_t diag_diff = std::abs(gap_s - gap_q);

                if (diag_diff <= max_gap) {
                    if (dpp[j] + 1 > dpp[i]) {
                        dpp[i] = dpp[j] + 1;
                        prevp[i] = static_cast<int32_t>(j);
                    }
                }
            }
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
            idx = (prev[idx] >= 0) ? static_cast<size_t>(prev[idx]) : SIZE_MAX;
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
