#include "search/result_dedup.hpp"

#include "index/ksx_reader.hpp"
#include "io/result_writer.hpp"
#include "search/parallel_search.hpp"

#include <algorithm>
#include <cstdint>

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
    const uint32_t shift      = ksx->fragment_start(h.cr.seq_id) - 1u;
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

}  // namespace ikafssn
