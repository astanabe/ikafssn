#pragma once

// Shared candidate-count limits.  Stage 1, Stage 2 and Stage 3 all cap their
// result sets the same way — a per-group top-N and a per-query in-total top-L —
// so the ranking rules live here once and each stage supplies its own group
// key and score.

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <vector>

namespace ikafssn {

struct ExtJob;
struct JobState;

// Keep at most `n` elements per group, ranked by score, in the manner of the
// -*_max_nhit_per_subject_mode switches:
//   tie_inclusive == false  keep the n highest, ties broken by arrival order
//   tie_inclusive == true   also keep every element tying the n-th score
// `group_less` / `group_eq` compare element indices by group key only; the
// ranking below appends the score (descending) and the arrival index
// (ascending).  `keep` is filled in input order.  n == 0 or cnt < 2 keeps
// everything.
template <typename GroupLess, typename GroupEq, typename ScoreOf>
void group_topn_keep(std::size_t cnt, uint32_t n, bool tie_inclusive,
                     GroupLess group_less, GroupEq group_eq, ScoreOf score_of,
                     std::vector<char>& keep) {
    if (n == 0 || cnt < 2) {
        keep.assign(cnt, 1);
        return;
    }
    keep.assign(cnt, 0);

    // One sort by (group, score desc, arrival asc) leaves each group as a
    // contiguous run that is already ranked.
    std::vector<uint32_t> order(cnt);
    for (std::size_t i = 0; i < cnt; ++i) order[i] = static_cast<uint32_t>(i);
    std::sort(order.begin(), order.end(), [&](uint32_t a, uint32_t b) {
        if (group_less(a, b)) return true;
        if (group_less(b, a)) return false;
        const auto sa = score_of(a);
        const auto sb = score_of(b);
        if (sa != sb) return sa > sb;
        return a < b;
    });

    const std::size_t nn = static_cast<std::size_t>(n);
    std::size_t i = 0;
    while (i < cnt) {
        std::size_t j = i + 1;
        while (j < cnt && group_eq(order[i], order[j])) ++j;
        const std::size_t run = j - i;
        if (run <= nn) {
            for (std::size_t t = i; t < j; ++t) keep[order[t]] = 1;
        } else if (!tie_inclusive) {
            for (std::size_t t = i; t < i + nn; ++t) keep[order[t]] = 1;
        } else {
            const auto s_n = score_of(order[i + nn - 1]);
            for (std::size_t t = i; t < j; ++t) {
                if (score_of(order[t]) < s_n) break;  // ranked descending
                keep[order[t]] = 1;
            }
        }
        i = j;
    }
}

// L-th highest score in `scores` (reorders it).  Returns the score type's
// lowest value when there are fewer than L, so that comparing `score >=
// threshold` keeps every element.
template <typename Score>
Score in_total_threshold(std::vector<Score>& scores, uint32_t L) {
    if (L == 0 || scores.size() < L) return std::numeric_limits<Score>::lowest();
    std::nth_element(scores.begin(), scores.begin() + (L - 1), scores.end(),
                     std::greater<Score>());
    return scores[L - 1];
}

// Per-group in_total_threshold.  `groups` and `scores` are parallel arrays and
// every group id is below `ngroups`.
template <typename Score>
std::vector<Score> in_total_thresholds(const std::vector<uint32_t>& groups,
                                       const std::vector<Score>& scores,
                                       std::size_t ngroups, uint32_t L) {
    std::vector<Score> thr(ngroups, std::numeric_limits<Score>::lowest());
    if (L == 0 || groups.empty()) return thr;

    // Bucket the scores by group into one flat buffer, so each group's run can
    // be selected in place.
    std::vector<uint32_t> begin(ngroups + 1, 0);
    for (uint32_t g : groups) ++begin[g + 1];
    for (std::size_t g = 0; g < ngroups; ++g) begin[g + 1] += begin[g];
    std::vector<Score> flat(groups.size());
    {
        std::vector<uint32_t> fill(begin.begin(), begin.end() - 1);
        for (std::size_t i = 0; i < groups.size(); ++i)
            flat[fill[groups[i]]++] = scores[i];
    }

    for (std::size_t g = 0; g < ngroups; ++g) {
        const std::size_t cnt = begin[g + 1] - begin[g];
        if (cnt < L) continue;  // fewer than L: every score is kept
        auto run = flat.begin() + begin[g];
        std::nth_element(run, run + (L - 1), run + cnt, std::greater<Score>());
        thr[g] = *(run + (L - 1));
    }
    return thr;
}

// Stage 1 in-total (L) limit for modes 2/3: keep, per query (across all volumes
// and strands), the top-L candidates plus every candidate tying the L-th
// coverscore, by pruning each ext_job's candidate list.  The prune is
// order-preserving, so the candidates stay ordered by ascending SeqId as
// Stage 2A requires.  L == 0 is a no-op.
void limit_candidates_in_total(std::vector<JobState>& states,
                               const std::vector<ExtJob>& ext_jobs, uint32_t L);

}  // namespace ikafssn
