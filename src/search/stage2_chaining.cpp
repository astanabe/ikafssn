#include "search/stage2_chaining.hpp"
#include "search/diagonal_filter.hpp"

#include <algorithm>
#include <cstdlib>

namespace ikafssn {

namespace {

// Per-thread working buffers for chain_hits().  They only ever grow, so once a
// thread has seen its largest subject the steady state allocates nothing.
struct ChainScratch {
    std::vector<Hit>      work;           // hits still available for extraction
    std::vector<Hit>      merged;         // natural-merge destination
    std::vector<uint32_t> q;              // DP input, split out of `work`
    std::vector<uint32_t> s;
    std::vector<int32_t>  diag;           // s - q, per hit
    std::vector<uint32_t> dp;             // longest chain ending at each hit
    std::vector<int32_t>  prev;           // that chain's predecessor, -1 if none
    std::vector<size_t>   chain_indices;  // traceback, strictly decreasing
    std::vector<int32_t>  filter_diags;   // diagonal_filter's working buffer
};

inline ChainScratch& tls_chain_scratch() {
    thread_local ChainScratch scratch;
    return scratch;
}

// Grow `v` to at least n elements and return its data.  Callers overwrite
// every element they read, so the grown part needs no particular value.
template <typename T>
inline T* scratch_ptr(std::vector<T>& v, size_t n) {
    if (v.size() < n) v.resize(n);
    return v.data();
}

// Function objects, not free functions: std::sort / std::unique would
// otherwise call the comparison through a function pointer.
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

    // Order by (q_pos, s_pos).  Stage 2A emits one monotonic run per seed
    // template, ascending on the forward strand and descending on the
    // reverse-complement strand (fwd_pos = len - pos - span), so detect and
    // merge the runs.  Reversing a run only reorders hits the dedup collapses.
    work.assign(raw_hits.begin(), raw_hits.end());
    {
        const size_t n = work.size();
        size_t i = 0, nruns = 0, run1_end = 0;
        while (i < n) {
            const size_t start = i++;
            if (i < n && HitLess{}(work[i], work[i - 1])) {
                while (i < n && HitLess{}(work[i], work[i - 1])) i++;
                std::reverse(work.begin() + start, work.begin() + i);
            } else {
                while (i < n && !HitLess{}(work[i], work[i - 1])) i++;
            }
            if (++nruns == 1) run1_end = i;
            if (nruns > 2) break;  // hand the rest to std::sort
        }
        if (nruns == 2) {
            Hit* const dst = scratch_ptr(scratch.merged, n);
            std::merge(work.begin(), work.begin() + run1_end,
                       work.begin() + run1_end, work.begin() + n,
                       dst, HitLess{});
            work.swap(scratch.merged);
            work.resize(n);
        } else if (nruns > 2) {
            std::sort(work.begin(), work.end(), HitLess{});
        }
    }

    // Drop the duplicate (q_pos, s_pos) pairs degenerate base expansion
    // produces; they would otherwise inflate the chainscore.
    work.erase(std::unique(work.begin(), work.end(), HitSamePos{}), work.end());

    diagonal_filter(work, config.min_nhit_diag, scratch.filter_diags);
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

        // Per-hit diagonal, precomputed: the gap test between hits i and j is
        //   |(s_i - s_j) - (q_i - q_j)| = |diag_i - diag_j|.
        // The buffers live in thread-local storage, so index them through
        // local pointers to keep the DP loop free of base-pointer reloads.
        const Hit* const src = work.data();
        uint32_t* const q = scratch_ptr(scratch.q, n);
        uint32_t* const s = scratch_ptr(scratch.s, n);
        int32_t* const diag = scratch_ptr(scratch.diag, n);
        uint32_t* const dp = scratch_ptr(scratch.dp, n);
        int32_t* const prev = scratch_ptr(scratch.prev, n);
        for (size_t i = 0; i < n; i++) {
            q[i] = src[i].q_pos;
            s[i] = src[i].s_pos;
            diag[i] = static_cast<int32_t>(src[i].s_pos) -
                      static_cast<int32_t>(src[i].q_pos);
        }

        const int64_t max_gap = static_cast<int64_t>(config.max_gap);
        const size_t lookback = config.chain_max_lookback;

        dp[0] = 1;
        prev[0] = -1;
        for (size_t i = 1; i < n; i++) {
            const size_t j_start = (lookback > 0 && i > lookback) ? (i - lookback) : 0;
            const uint32_t q_i = q[i];
            const uint32_t s_i = s[i];
            const int64_t diag_i = diag[i];
            // Ascending j, and only a strictly longer chain replaces the
            // incumbent, so a tie keeps the smallest j.  Changing the scan
            // order would change the reported coordinates.
            uint32_t dp_i = 1;
            int32_t prev_i = -1;
            for (size_t j = j_start; j < i; j++) {
                if (q[j] >= q_i) continue;
                if (s[j] >= s_i) continue;
                if (dp[j] < dp_i) continue;
                if (std::llabs(diag_i - static_cast<int64_t>(diag[j])) > max_gap) continue;
                dp_i = dp[j] + 1;
                prev_i = static_cast<int32_t>(j);
            }
            dp[i] = dp_i;
            prev[i] = prev_i;
        }

        // Find best chain endpoint
        size_t best_idx = 0;
        for (size_t i = 1; i < n; i++) {
            if (dp[i] > dp[best_idx]) {
                best_idx = i;
            }
        }

        uint32_t best_score = dp[best_idx];
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

        // Strict take-1: no next iteration, so skip the removal pass.
        if (!unlimited && !tie_inclusive && n_limit == 1) break;

        // Drop the chain's hits.  prev[i] < i for every link, so
        // chain_indices is strictly decreasing and walking it backwards gives
        // the removal targets in ascending order — one compaction pass.
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
