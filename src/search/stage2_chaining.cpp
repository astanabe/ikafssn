#include "search/stage2_chaining.hpp"
#include "search/diagonal_filter.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>

namespace ikafssn {

ChainResult chain_hits(const std::vector<Hit>& raw_hits,
                       SeqId seq_id,
                       int k,
                       bool is_reverse,
                       const Stage2Config& config) {
    ChainResult result{};
    result.seq_id = seq_id;
    result.score = 0;
    result.is_reverse = is_reverse;

    if (raw_hits.empty()) return result;

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
    if (hits.empty()) return result;

    // Already sorted by (q_pos, s_pos) from dedup step

    // Step 3: O(n^2) chaining DP
    size_t n = hits.size();
    std::vector<uint32_t> dp(n, 1);    // each hit is a chain of length 1
    std::vector<int32_t> prev(n, -1);  // traceback pointers

    for (size_t i = 1; i < n; i++) {
        size_t j_start = (config.chain_max_lookback > 0 && i > config.chain_max_lookback)
                          ? (i - config.chain_max_lookback) : 0;
        for (size_t j = j_start; j < i; j++) {
            // Collinearity constraint: q_pos[j] < q_pos[i] and s_pos[j] < s_pos[i]
            // Both must strictly increase to ensure each chain element comes
            // from a distinct query k-mer position.
            if (hits[j].q_pos >= hits[i].q_pos) continue;
            if (hits[j].s_pos >= hits[i].s_pos) continue;

            // Gap constraint
            int64_t gap_q = static_cast<int64_t>(hits[i].q_pos) - static_cast<int64_t>(hits[j].q_pos);
            int64_t gap_s = static_cast<int64_t>(hits[i].s_pos) - static_cast<int64_t>(hits[j].s_pos);
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
    if (best_score < config.min_score) return result;

    // Step 4: traceback to find chain start
    size_t chain_start_idx = best_idx;
    while (prev[chain_start_idx] >= 0) {
        chain_start_idx = static_cast<size_t>(prev[chain_start_idx]);
    }

    result.score = best_score;
    result.q_start = hits[chain_start_idx].q_pos;
    result.q_end = hits[best_idx].q_pos + static_cast<uint32_t>(k);
    result.s_start = hits[chain_start_idx].s_pos;
    result.s_end = hits[best_idx].s_pos + static_cast<uint32_t>(k);

    return result;
}

} // namespace ikafssn
