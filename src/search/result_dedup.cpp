#include "search/result_dedup.hpp"

#include "index/ksx_reader.hpp"
#include "io/result_writer.hpp"
#include "search/hit_limits.hpp"
#include "search/parallel_search.hpp"
#include "search/parent_hit.hpp"
#include "util/query_order.hpp"

#include <algorithm>
#include <climits>
#include <cstdint>
#include <string>
#include <unordered_map>

namespace ikafssn {

namespace {

struct StageTwoKey {
    size_t   query_idx;
    size_t   volume_idx;
    uint32_t parent_idx;
    uint32_t parent_sstart;
    uint32_t parent_send;
    uint32_t chainscore;
    bool     is_reverse;
};

inline StageTwoKey make_stage2_key(
    const OrchestratorHit& h,
    const std::vector<const KsxReader*>& ksx_per_volume) {
    const auto* ksx = ksx_per_volume[h.volume_idx];
    const uint32_t parent_idx = ksx->parent_index(h.cr.seq_id);
    const uint32_t shift      = fragment_shift(*ksx, h.cr.seq_id);
    return StageTwoKey{
        h.query_idx,
        h.volume_idx,
        parent_idx,
        h.cr.s_start + shift,
        h.cr.s_end   + shift,
        h.cr.chainscore,
        h.cr.is_reverse,
    };
}

inline bool stage2_key_less(const StageTwoKey& a, const StageTwoKey& b) {
    if (a.query_idx != b.query_idx)         return a.query_idx     < b.query_idx;
    if (a.volume_idx != b.volume_idx)       return a.volume_idx    < b.volume_idx;
    if (a.parent_idx != b.parent_idx)       return a.parent_idx    < b.parent_idx;
    if (a.is_reverse != b.is_reverse)       return a.is_reverse    < b.is_reverse;
    if (a.parent_sstart != b.parent_sstart) return a.parent_sstart < b.parent_sstart;
    if (a.parent_send != b.parent_send)     return a.parent_send   < b.parent_send;
    return a.chainscore < b.chainscore;
}

inline bool stage2_key_eq(const StageTwoKey& a, const StageTwoKey& b) {
    return a.query_idx     == b.query_idx
        && a.volume_idx    == b.volume_idx
        && a.parent_idx    == b.parent_idx
        && a.is_reverse    == b.is_reverse
        && a.parent_sstart == b.parent_sstart
        && a.parent_send   == b.parent_send
        && a.chainscore    == b.chainscore;
}

}  // namespace

void dedup_stage2_orchestrator_hits(
    std::vector<OrchestratorHit>& hits,
    const std::vector<const KsxReader*>& ksx_per_volume) {
    if (hits.size() < 2) return;

    const size_t n = hits.size();
    std::vector<StageTwoKey> keys(n);
    std::vector<size_t> idx(n);
    for (size_t i = 0; i < n; ++i) {
        keys[i] = make_stage2_key(hits[i], ksx_per_volume);
        idx[i]  = i;
    }

    std::sort(idx.begin(), idx.end(),
              [&](size_t a, size_t b) {
                  return stage2_key_less(keys[a], keys[b]);
              });

    std::vector<OrchestratorHit> out;
    out.reserve(n);
    const StageTwoKey* prev = nullptr;
    for (size_t i : idx) {
        if (prev != nullptr && stage2_key_eq(*prev, keys[i])) continue;
        out.push_back(std::move(hits[i]));
        prev = &keys[i];
    }
    hits = std::move(out);
}

void select_parent_topn(
    std::vector<OrchestratorHit>& hits,
    uint32_t n,
    uint8_t mode,
    const std::vector<const KsxReader*>& ksx_per_volume) {
    if (n == 0 || hits.size() < 2) return;

    const bool strand_split   = (mode == 2 || mode == 4);
    const bool tie_inclusive  = (mode == 3 || mode == 4);

    struct ParentKey {
        size_t   query_idx;
        size_t   volume_idx;
        uint32_t parent_idx;
        bool     is_reverse;  // meaningful only when strand_split
    };

    const size_t cnt = hits.size();
    std::vector<ParentKey> keys(cnt);
    for (size_t i = 0; i < cnt; ++i) {
        const auto& h = hits[i];
        const auto* ksx = ksx_per_volume[h.volume_idx];
        keys[i] = ParentKey{
            h.query_idx,
            h.volume_idx,
            ksx->parent_index(h.cr.seq_id),
            strand_split ? h.cr.is_reverse : false,
        };
    }

    std::vector<char> keep;
    group_topn_keep(
        cnt, n, tie_inclusive,
        [&](uint32_t a, uint32_t b) {
            const ParentKey& ka = keys[a];
            const ParentKey& kb = keys[b];
            if (ka.query_idx  != kb.query_idx)  return ka.query_idx  < kb.query_idx;
            if (ka.volume_idx != kb.volume_idx) return ka.volume_idx < kb.volume_idx;
            if (ka.parent_idx != kb.parent_idx) return ka.parent_idx < kb.parent_idx;
            return ka.is_reverse < kb.is_reverse;
        },
        [&](uint32_t a, uint32_t b) {
            const ParentKey& ka = keys[a];
            const ParentKey& kb = keys[b];
            return ka.query_idx  == kb.query_idx
                && ka.volume_idx == kb.volume_idx
                && ka.parent_idx == kb.parent_idx
                && ka.is_reverse == kb.is_reverse;
        },
        [&](uint32_t a) { return hits[a].cr.chainscore; },
        keep);

    std::vector<OrchestratorHit> out;
    out.reserve(cnt);
    for (size_t i = 0; i < cnt; ++i)
        if (keep[i]) out.push_back(std::move(hits[i]));
    hits = std::move(out);
}

void dedup_stage3_output_hits(std::vector<OutputHit>& hits) {
    if (hits.size() < 2) return;

    // Skip-marker rows are not deduplicated.  Pull them aside so the dedup
    // pass runs only over real hits, and re-append them at the end in their
    // original order.
    std::vector<OutputHit> skip;
    std::vector<OutputHit> real;
    skip.reserve(hits.size());
    real.reserve(hits.size());
    for (auto& h : hits) {
        if (h.skip_reason != 0) skip.push_back(std::move(h));
        else                    real.push_back(std::move(h));
    }

    std::sort(real.begin(), real.end(),
              [](const OutputHit& a, const OutputHit& b) {
                  if (a.qseqid != b.qseqid)   return a.qseqid  < b.qseqid;
                  if (a.sseqid != b.sseqid)   return a.sseqid  < b.sseqid;
                  if (a.sstrand != b.sstrand) return a.sstrand < b.sstrand;
                  if (a.send != b.send)       return a.send    < b.send;
                  return a.alnscore < b.alnscore;
              });
    auto last = std::unique(real.begin(), real.end(),
              [](const OutputHit& a, const OutputHit& b) {
                  return a.qseqid  == b.qseqid
                      && a.sseqid  == b.sseqid
                      && a.sstrand == b.sstrand
                      && a.send    == b.send
                      && a.alnscore == b.alnscore;
              });
    real.erase(last, real.end());

    hits.clear();
    hits.reserve(real.size() + skip.size());
    for (auto& h : real) hits.push_back(std::move(h));
    for (auto& h : skip) hits.push_back(std::move(h));
}

void select_parent_topn_output(std::vector<OutputHit>& hits,
                               uint32_t n,
                               uint8_t mode) {
    if (n == 0 || hits.size() < 2) return;

    const bool strand_split  = (mode == 2 || mode == 4);
    const bool tie_inclusive = (mode == 3 || mode == 4);

    const size_t cnt = hits.size();
    std::vector<char> keep;
    group_topn_keep(
        cnt, n, tie_inclusive,
        [&](uint32_t a, uint32_t b) {
            const OutputHit& ha = hits[a];
            const OutputHit& hb = hits[b];
            if (ha.qseqid != hb.qseqid) return ha.qseqid < hb.qseqid;
            if (ha.sseqid != hb.sseqid) return ha.sseqid < hb.sseqid;
            if (strand_split) return ha.sstrand < hb.sstrand;
            return false;
        },
        [&](uint32_t a, uint32_t b) {
            const OutputHit& ha = hits[a];
            const OutputHit& hb = hits[b];
            return ha.qseqid == hb.qseqid
                && ha.sseqid == hb.sseqid
                && (!strand_split || ha.sstrand == hb.sstrand);
        },
        [&](uint32_t a) { return hits[a].alnscore; },
        keep);

    std::vector<OutputHit> out;
    out.reserve(cnt);
    for (size_t i = 0; i < cnt; ++i)
        if (keep[i]) out.push_back(std::move(hits[i]));
    hits = std::move(out);
}

void apply_in_total_output(std::vector<OutputHit>& hits, uint32_t L) {
    if (L == 0 || hits.empty()) return;

    // Per-qseqid threshold = the L-th highest alnscore, then keep every hit at
    // or above it.  A query with at most L hits gets INT32_MIN and keeps all.
    QueryOrder qorder;
    std::vector<uint32_t> groups(hits.size());
    std::vector<int32_t> scores(hits.size());
    for (size_t i = 0; i < hits.size(); ++i) {
        groups[i] = qorder.id_of(hits[i].qseqid);
        scores[i] = hits[i].alnscore;
    }
    const auto thr = in_total_thresholds(groups, scores,
                                         qorder.qseqids().size(), L);

    std::vector<OutputHit> out;
    out.reserve(hits.size());
    for (size_t i = 0; i < hits.size(); ++i) {
        if (hits[i].alnscore >= thr[groups[i]]) out.push_back(std::move(hits[i]));
    }
    hits = std::move(out);
}

}  // namespace ikafssn
