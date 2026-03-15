#include "search/stage2_chaining.hpp"
#include "search/diagonal_filter.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>

namespace ikafssn {

std::vector<ChainResult> chain_hits(const std::vector<Hit>& raw_hits,
                                    SeqId seq_id,
                                    int k,
                                    bool is_reverse,
                                    const Stage2Config& config) {
    if (raw_hits.empty()) return {};

    // Step 1: deduplicate (q_pos, s_pos) pairs from degenerate base expansion
    std::vector<Hit> deduped = raw_hits;
    std::sort(deduped.begin(), deduped.end(), [](const Hit& a, const Hit& b) {
        return a.q_pos < b.q_pos || (a.q_pos == b.q_pos && a.s_pos < b.s_pos);
    });
    deduped.erase(std::unique(deduped.begin(), deduped.end(),
        [](const Hit& a, const Hit& b) {
            return a.q_pos == b.q_pos && a.s_pos == b.s_pos;
        }), deduped.end());

    // Step 2: diagonal filter
    std::vector<Hit> hits = diagonal_filter(deduped, config.min_diag_hits);
    if (hits.empty()) return {};

    // Already sorted by (q_pos, s_pos) from dedup step

    // Determine max iterations
    uint32_t max_chains = config.max_nhit_per_subject;
    if (max_chains == 0) max_chains = UINT32_MAX; // unlimited

    std::vector<ChainResult> results;

    // Working copy of hits for iterative removal
    std::vector<Hit> remaining = std::move(hits);

    for (uint32_t iter = 0; iter < max_chains; iter++) {
        size_t n = remaining.size();
        if (n == 0) break;

        // DP arrays
        std::vector<uint32_t> dp(n, 1);
        std::vector<int32_t> prev(n, -1);

        for (size_t i = 1; i < n; i++) {
            size_t j_start = (config.chain_max_lookback > 0 && i > config.chain_max_lookback)
                              ? (i - config.chain_max_lookback) : 0;
            for (size_t j = j_start; j < i; j++) {
                if (remaining[j].q_pos >= remaining[i].q_pos) continue;
                if (remaining[j].s_pos >= remaining[i].s_pos) continue;

                int64_t gap_q = static_cast<int64_t>(remaining[i].q_pos) - static_cast<int64_t>(remaining[j].q_pos);
                int64_t gap_s = static_cast<int64_t>(remaining[i].s_pos) - static_cast<int64_t>(remaining[j].s_pos);
                int64_t diag_diff = std::abs(gap_s - gap_q);

                if (diag_diff <= static_cast<int64_t>(config.max_gap)) {
                    if (dp[j] + 1 > dp[i]) {
                        dp[i] = dp[j] + 1;
                        prev[i] = static_cast<int32_t>(j);
                    }
                }
            }
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

        // Traceback to collect chain hit indices
        std::vector<size_t> chain_indices;
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
        cr.q_start = remaining[chain_start_idx].q_pos;
        cr.q_end = remaining[chain_end_idx].q_pos + static_cast<uint32_t>(k);
        cr.s_start = remaining[chain_start_idx].s_pos;
        cr.s_end = remaining[chain_end_idx].s_pos + static_cast<uint32_t>(k);

        results.push_back(cr);

        // Early return for single chain (default path, no removal overhead)
        if (max_chains == 1) break;

        // Remove used hits for next iteration
        std::vector<bool> used(n, false);
        for (size_t ci : chain_indices) {
            used[ci] = true;
        }

        std::vector<Hit> next_remaining;
        next_remaining.reserve(n - chain_indices.size());
        for (size_t i = 0; i < n; i++) {
            if (!used[i]) {
                next_remaining.push_back(remaining[i]);
            }
        }
        remaining = std::move(next_remaining);
    }

    return results;
}

} // namespace ikafssn
