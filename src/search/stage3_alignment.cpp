#include "search/stage3_alignment.hpp"
#include "core/spaced_seed.hpp"
#include "io/blastdb_reader.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include <parasail.h>

namespace ikafssn {

// Walk a parasail CIGAR to compute npositive, nnegative, and build a CIGAR string.
struct CigarStats {
    uint32_t npositive = 0;
    uint32_t nnegative = 0;
    uint32_t aln_len = 0;
    std::string cigar_str;
};

static CigarStats walk_cigar(const parasail_cigar_t* cigar) {
    CigarStats stats;
    for (int i = 0; i < cigar->len; i++) {
        char op = parasail_cigar_decode_op(cigar->seq[i]);
        uint32_t len = parasail_cigar_decode_len(cigar->seq[i]);
        stats.aln_len += len;
        switch (op) {
            case '=': stats.npositive += len; break;
            case 'X': stats.nnegative += len; break;
            case 'I': break; // insertion in query
            case 'D': break; // deletion in query (insertion in ref)
            default: break;
        }
        stats.cigar_str += std::to_string(len);
        stats.cigar_str += op;
    }
    return stats;
}

// ---------------------------------------------------------------------------
// Stage 3 batching
//
// Hits sharing the same (qseqid, sseqid, sstrand) form an atomic group: the
// overlap-resolution loop assumes all such hits are visible at once, so a
// group cannot be split across batches.  Groups are bin-packed greedily
// against `posting_budget` (the residual heap budget already plumbed from
// Stage 1/2 madvise accounting).  A group whose own cost exceeds the budget
// is processed in a solo batch (over-budget).
// ---------------------------------------------------------------------------

static uint64_t compute_hit_cost(const OutputHit& h,
                                 const Stage3Config& cfg,
                                 bool ctx_is_ratio,
                                 double ctx_ratio,
                                 uint32_t ctx_abs)
{
    uint32_t ctx = ctx_is_ratio
        ? static_cast<uint32_t>(h.qlen * ctx_ratio)
        : ctx_abs;
    uint64_t subseq = static_cast<uint64_t>(h.send - h.sstart + 1) + 2ull * ctx;
    uint64_t tb = cfg.traceback
        ? 2ull * std::max<uint64_t>(h.qlen, subseq)
        : 0ull;
    return subseq + tb + 256; // const overhead: OutputHit + std::string headers + pad
}

struct Stage3Group {
    std::vector<size_t> indices;
    uint64_t cost = 0;
};

struct Stage3Batch {
    std::vector<size_t> group_idxs; // indices into the ordered groups vector
    uint64_t cost = 0;
    bool solo_oversize = false;
};

static std::vector<Stage3Batch>
plan_stage3_batches(const std::vector<Stage3Group>& groups_ordered,
                    uint64_t budget)
{
    std::vector<Stage3Batch> batches;
    Stage3Batch cur;
    auto flush = [&]() {
        if (!cur.group_idxs.empty()) {
            batches.push_back(std::move(cur));
            cur = Stage3Batch{};
        }
    };
    for (size_t gi = 0; gi < groups_ordered.size(); gi++) {
        const auto& g = groups_ordered[gi];
        if (budget > 0 && g.cost > budget) {
            flush();
            Stage3Batch solo;
            solo.group_idxs.push_back(gi);
            solo.cost = g.cost;
            solo.solo_oversize = true;
            batches.push_back(std::move(solo));
            continue;
        }
        if (budget > 0 && cur.cost + g.cost > budget && !cur.group_idxs.empty()) {
            flush();
        }
        cur.group_idxs.push_back(gi);
        cur.cost += g.cost;
    }
    flush();
    return batches;
}

std::vector<OutputHit> run_stage3(
    std::vector<OutputHit>& hits,
    const std::vector<FastaRecord>& queries,
    const std::string& db_path,
    const Stage3Config& config,
    bool context_is_ratio,
    double context_ratio,
    uint32_t context_abs,
    const Logger& logger)
{
    if (hits.empty()) return {};

    // Open BLAST DB volumes
    auto vol_paths = BlastDbReader::find_volume_paths(db_path);
    if (vol_paths.empty()) {
        logger.error("Stage 3: no BLAST DB volumes found at '%s'", db_path.c_str());
        return {};
    }

    std::vector<BlastDbReader> readers(vol_paths.size());
    for (size_t vi = 0; vi < vol_paths.size(); vi++) {
        if (!readers[vi].open(vol_paths[vi])) {
            logger.error("Stage 3: cannot open volume '%s'", vol_paths[vi].c_str());
            return {};
        }
    }

    // Build query lookup: query_id -> index in queries[]
    std::unordered_map<std::string, size_t> query_map;
    for (size_t i = 0; i < queries.size(); i++) {
        query_map[queries[i].id] = i;
    }

    // Validate hit volumes
    std::vector<bool> hit_valid(hits.size(), true);
    for (size_t i = 0; i < hits.size(); i++) {
        if (hits[i].volume >= readers.size()) {
            logger.warn("Stage 3: hit volume %u out of range (max %zu), skipping",
                        static_cast<unsigned>(hits[i].volume), readers.size());
            hit_valid[i] = false;
        }
    }

    // Build profiles up front (per unique query_idx x strand).  Profiles
    // are small and shared across all batches.
    const parasail_matrix_t* matrix = parasail_matrix_lookup(config.score_matrix.c_str());
    if (!matrix) {
        logger.error("Stage 3: unknown score matrix '%s'", config.score_matrix.c_str());
        return {};
    }

    struct ProfileEntry {
        parasail_profile_t* profile = nullptr;
        std::string seq; // keep alive for profile lifetime
    };
    std::unordered_map<std::string, ProfileEntry> profiles; // key: "qidx:strand"

    for (size_t i = 0; i < hits.size(); i++) {
        if (!hit_valid[i]) continue;
        auto qit = query_map.find(hits[i].qseqid);
        if (qit == query_map.end()) continue;
        size_t qi = qit->second;
        bool is_rev = (hits[i].sstrand == '-');
        std::string key = std::to_string(qi) + ":" + (is_rev ? "1" : "0");
        if (profiles.find(key) == profiles.end()) {
            ProfileEntry pe;
            pe.seq = is_rev ? reverse_complement_string(queries[qi].sequence)
                            : queries[qi].sequence;
            pe.profile = parasail_profile_create_sat(
                pe.seq.c_str(), static_cast<int>(pe.seq.size()), matrix);
            profiles[key] = std::move(pe);
        }
    }

    // Build atomic groups keyed by (qseqid, sseqid, sstrand).
    // The overlap-resolution loop assumes each group is visible whole, so
    // groups are the unit of batching.  Groups are stored in a parallel
    // vector so plan_stage3_batches can index them.
    std::unordered_map<std::string, size_t> group_index;
    std::vector<Stage3Group> groups;
    groups.reserve(64);
    for (size_t i = 0; i < hits.size(); i++) {
        if (!hit_valid[i]) continue;
        if (query_map.find(hits[i].qseqid) == query_map.end()) {
            // Skip hits whose query is missing — they cannot be aligned.
            hit_valid[i] = false;
            continue;
        }
        std::string key = hits[i].qseqid + "\t" + hits[i].sseqid + "\t" + hits[i].sstrand;
        auto it = group_index.find(key);
        if (it == group_index.end()) {
            group_index[key] = groups.size();
            groups.push_back(Stage3Group{});
            groups.back().indices.push_back(i);
        } else {
            groups[it->second].indices.push_back(i);
        }
    }

    for (auto& g : groups) {
        // Sort by sstart for the overlap-resolution loop's invariant.
        std::sort(g.indices.begin(), g.indices.end(),
            [&hits](size_t a, size_t b) { return hits[a].sstart < hits[b].sstart; });
        for (size_t idx : g.indices) {
            g.cost += compute_hit_cost(hits[idx], config,
                                       context_is_ratio, context_ratio, context_abs);
        }
    }

    // Plan batches.
    auto batches = plan_stage3_batches(groups, config.posting_budget);

    {
        size_t n_solo = 0;
        for (const auto& b : batches) if (b.solo_oversize) ++n_solo;
        logger.info("Stage 3 batch plan: %zu batch(es) over %zu group(s) "
                    "(budget=%llu, solo_oversize=%zu)",
                    batches.size(), groups.size(),
                    static_cast<unsigned long long>(config.posting_budget),
                    n_solo);
    }

    // Per-hit scratch storage.  Sized to hits.size() so we can index by
    // the original hit index, but each batch only fills entries for its
    // own hits and clears them when the batch ends.
    std::vector<std::string> subject_subseqs(hits.size());
    std::vector<uint32_t> ext_starts(hits.size(), 0);

    std::vector<OutputHit> filtered;
    filtered.reserve(hits.size());

    bool has_context = context_is_ratio ? (context_ratio > 0) : (context_abs > 0);

    // Per-batch loop.
    for (size_t bi = 0; bi < batches.size(); bi++) {
        const auto& batch = batches[bi];

        // Collect this batch's hit indices grouped by volume for
        // sequential mmap access.
        std::vector<std::vector<size_t>> hits_by_reader(readers.size());
        size_t batch_hit_count = 0;
        for (size_t gi : batch.group_idxs) {
            for (size_t hidx : groups[gi].indices) {
                if (!hit_valid[hidx]) continue;
                hits_by_reader[hits[hidx].volume].push_back(hidx);
                ++batch_hit_count;
            }
        }
        for (auto& vh : hits_by_reader) {
            std::sort(vh.begin(), vh.end(),
                [&hits](size_t a, size_t b) { return hits[a].oid < hits[b].oid; });
        }

        logger.debug("Stage 3 batch %zu/%zu: %zu groups, %zu hits, est=%llu bytes%s",
                     bi + 1, batches.size(),
                     batch.group_idxs.size(), batch_hit_count,
                     static_cast<unsigned long long>(batch.cost),
                     batch.solo_oversize ? " [solo oversize]" : "");
        if (batch.solo_oversize) {
            logger.warn("Stage 3 batch %zu/%zu: single group exceeds posting_budget "
                        "(%llu > %llu); processing solo (over-budget)",
                        bi + 1, batches.size(),
                        static_cast<unsigned long long>(batch.cost),
                        static_cast<unsigned long long>(config.posting_budget));
        }

        // Fetch subseqs.  Flatten the per-volume OID-sorted hit lists
        // into one ordered index sequence and walk it with a hit-parallel
        // parallel_for; iterations within a task stay in (volume, OID) order
        // for sequential mmap locality.  get_subsequence is lock-free, so
        // concurrent calls are safe.
        std::vector<size_t> ordered_hits;
        ordered_hits.reserve(batch_hit_count);
        for (size_t ri = 0; ri < readers.size(); ri++) {
            for (size_t h : hits_by_reader[ri]) ordered_hits.push_back(h);
        }
        for (size_t ri = 0; ri < readers.size(); ri++) {
            if (hits_by_reader[ri].empty()) continue;
            readers[ri].set_mmap_strategy(BlastDbReader::MMapStrategy::kNormal);
        }
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, ordered_hits.size(), 16),
            [&](const tbb::blocked_range<size_t>& r) {
                for (size_t i = r.begin(); i < r.end(); i++) {
                    size_t hit_idx = ordered_hits[i];
                    uint16_t vol = hits[hit_idx].volume;
                    uint32_t oid = hits[hit_idx].oid;
                    uint32_t seq_len = readers[vol].seq_length(oid);

                    uint32_t query_len = 0;
                    auto qit = query_map.find(hits[hit_idx].qseqid);
                    if (qit != query_map.end()) {
                        query_len = static_cast<uint32_t>(queries[qit->second].sequence.size());
                    }

                    uint32_t ctx = context_is_ratio
                        ? static_cast<uint32_t>(query_len * context_ratio)
                        : context_abs;

                    uint32_t ext_start = (hits[hit_idx].sstart >= ctx)
                        ? hits[hit_idx].sstart - ctx : 0;
                    uint32_t ext_end = std::min(hits[hit_idx].send + ctx, seq_len - 1);

                    subject_subseqs[hit_idx] = readers[vol].get_subsequence(
                        oid, ext_start, ext_end);
                    ext_starts[hit_idx] = ext_start;
                    hits[hit_idx].slen = seq_len;
                }
            });
        // Subseqs are now in `subject_subseqs` (heap-owned strings); the
        // alignment step does not touch the mmap further this batch, so
        // release the cached pages before the next batch's fetch.
        for (size_t ri = 0; ri < readers.size(); ri++) {
            if (hits_by_reader[ri].empty()) continue;
            readers[ri].set_mmap_strategy(BlastDbReader::MMapStrategy::kDontNeed);
        }

        // Collect valid hit indices for parallel alignment.
        std::vector<size_t> valid_indices;
        valid_indices.reserve(batch_hit_count);
        for (size_t gi : batch.group_idxs) {
            for (size_t hidx : groups[gi].indices) {
                if (hit_valid[hidx] && !subject_subseqs[hidx].empty()) {
                    valid_indices.push_back(hidx);
                }
            }
        }

        // Parallel alignment.
        tbb::parallel_for(size_t(0), valid_indices.size(), [&](size_t vi) {
            size_t idx = valid_indices[vi];
            auto qit = query_map.find(hits[idx].qseqid);
            size_t qi = qit->second;
            bool is_rev = (hits[idx].sstrand == '-');
            std::string key = std::to_string(qi) + ":" + (is_rev ? "1" : "0");
            const auto& pe = profiles.at(key);

            const char* subj = subject_subseqs[idx].c_str();
            int slen = static_cast<int>(subject_subseqs[idx].size());

            if (config.traceback) {
                parasail_result_t* result = parasail_sg_trace_striped_profile_sat(
                    pe.profile, subj, slen, config.gapopen, config.gapext);

                hits[idx].alnscore = result->score;

                parasail_cigar_t* cigar = parasail_result_get_cigar(
                    result, pe.seq.c_str(), static_cast<int>(pe.seq.size()),
                    subj, slen, matrix);

                hits[idx].qstart = static_cast<uint32_t>(cigar->beg_query);
                hits[idx].qend = static_cast<uint32_t>(result->end_query);
                hits[idx].sstart = ext_starts[idx] + static_cast<uint32_t>(cigar->beg_ref);
                hits[idx].send = ext_starts[idx] + static_cast<uint32_t>(result->end_ref);

                CigarStats cs = walk_cigar(cigar);
                hits[idx].npositive = cs.npositive;
                hits[idx].nnegative = cs.nnegative;
                hits[idx].cigar = cs.cigar_str;
                hits[idx].ppositive = (cs.aln_len > 0) ? 100.0 * cs.npositive / cs.aln_len : 0.0;

                parasail_traceback_t* tb = parasail_result_get_traceback(
                    result, pe.seq.c_str(), static_cast<int>(pe.seq.size()),
                    subj, slen, matrix, '|', '*', ' ');
                if (tb) {
                    hits[idx].qseq = tb->query;
                    hits[idx].sseq = tb->ref;
                    parasail_traceback_free(tb);
                }

                parasail_cigar_free(cigar);
                parasail_result_free(result);
            } else {
                parasail_result_t* result = parasail_sg_striped_profile_sat(
                    pe.profile, subj, slen, config.gapopen, config.gapext);

                hits[idx].alnscore = result->score;
                hits[idx].qend = static_cast<uint32_t>(result->end_query);
                hits[idx].send = ext_starts[idx] + static_cast<uint32_t>(result->end_ref);

                parasail_result_free(result);
            }
        });

        // Overlap resolution for multi-chain hits in this batch's groups
        // (context > 0 only).  Groups in different batches share no
        // state, so the loop is naturally batch-local.
        if (has_context) {
            for (size_t gi : batch.group_idxs) {
                auto& group = groups[gi];
                if (group.indices.size() < 2) continue;

                bool changed = true;
                while (changed) {
                    changed = false;

                    for (size_t pi = 0; pi + 1 < group.indices.size(); pi++) {
                        size_t idx_a = group.indices[pi];
                        size_t idx_b = group.indices[pi + 1];

                        if (!hit_valid[idx_a] || !hit_valid[idx_b]) continue;
                        if (hits[idx_a].send < hits[idx_b].sstart) continue;

                        size_t keep_idx, clamp_idx;
                        if (hits[idx_a].chainscore >= hits[idx_b].chainscore) {
                            keep_idx = idx_a;
                            clamp_idx = idx_b;
                        } else {
                            keep_idx = idx_b;
                            clamp_idx = idx_a;
                        }

                        uint32_t new_ext_start, new_ext_end;
                        uint32_t oid = hits[clamp_idx].oid;
                        uint32_t seq_len = hits[clamp_idx].slen;
                        if (seq_len == 0) seq_len = 1;

                        auto qit2 = query_map.find(hits[clamp_idx].qseqid);
                        if (qit2 == query_map.end()) continue;
                        uint32_t query_len2 = static_cast<uint32_t>(queries[qit2->second].sequence.size());
                        uint32_t ctx2 = context_is_ratio
                            ? static_cast<uint32_t>(query_len2 * context_ratio)
                            : context_abs;

                        if (hits[clamp_idx].sstart <= hits[keep_idx].sstart) {
                            uint32_t boundary = hits[keep_idx].sstart;
                            if (hits[clamp_idx].sstart >= boundary) {
                                hit_valid[clamp_idx] = false;
                                changed = true;
                                continue;
                            }
                            new_ext_start = ext_starts[clamp_idx];
                            new_ext_end = (boundary > 0) ? boundary - 1 : 0;
                        } else {
                            uint32_t boundary = hits[keep_idx].send + 1;
                            if (boundary >= hits[clamp_idx].send) {
                                hit_valid[clamp_idx] = false;
                                changed = true;
                                continue;
                            }
                            new_ext_start = boundary;
                            new_ext_end = std::min(hits[clamp_idx].send + ctx2, seq_len - 1);
                        }

                        if (new_ext_start >= new_ext_end) {
                            hit_valid[clamp_idx] = false;
                            changed = true;
                            continue;
                        }

                        uint16_t vol = hits[clamp_idx].volume;
                        if (vol >= readers.size()) {
                            hit_valid[clamp_idx] = false;
                            changed = true;
                            continue;
                        }
                        subject_subseqs[clamp_idx] = readers[vol].get_subsequence(
                            oid, new_ext_start, new_ext_end);
                        ext_starts[clamp_idx] = new_ext_start;

                        size_t qi2 = qit2->second;
                        bool is_rev2 = (hits[clamp_idx].sstrand == '-');
                        std::string pkey = std::to_string(qi2) + ":" + (is_rev2 ? "1" : "0");
                        auto pit = profiles.find(pkey);
                        if (pit == profiles.end()) continue;
                        const auto& pe2 = pit->second;

                        const char* subj2 = subject_subseqs[clamp_idx].c_str();
                        int slen2 = static_cast<int>(subject_subseqs[clamp_idx].size());

                        if (config.traceback) {
                            parasail_result_t* result2 = parasail_sg_trace_striped_profile_sat(
                                pe2.profile, subj2, slen2, config.gapopen, config.gapext);
                            hits[clamp_idx].alnscore = result2->score;

                            parasail_cigar_t* cigar2 = parasail_result_get_cigar(
                                result2, pe2.seq.c_str(), static_cast<int>(pe2.seq.size()),
                                subj2, slen2, matrix);
                            hits[clamp_idx].qstart = static_cast<uint32_t>(cigar2->beg_query);
                            hits[clamp_idx].qend = static_cast<uint32_t>(result2->end_query);
                            hits[clamp_idx].sstart = new_ext_start + static_cast<uint32_t>(cigar2->beg_ref);
                            hits[clamp_idx].send = new_ext_start + static_cast<uint32_t>(result2->end_ref);

                            CigarStats cs2 = walk_cigar(cigar2);
                            hits[clamp_idx].npositive = cs2.npositive;
                            hits[clamp_idx].nnegative = cs2.nnegative;
                            hits[clamp_idx].cigar = cs2.cigar_str;
                            hits[clamp_idx].ppositive = (cs2.aln_len > 0) ? 100.0 * cs2.npositive / cs2.aln_len : 0.0;

                            parasail_traceback_t* tb2 = parasail_result_get_traceback(
                                result2, pe2.seq.c_str(), static_cast<int>(pe2.seq.size()),
                                subj2, slen2, matrix, '|', '*', ' ');
                            if (tb2) {
                                hits[clamp_idx].qseq = tb2->query;
                                hits[clamp_idx].sseq = tb2->ref;
                                parasail_traceback_free(tb2);
                            }

                            parasail_cigar_free(cigar2);
                            parasail_result_free(result2);
                        } else {
                            parasail_result_t* result2 = parasail_sg_striped_profile_sat(
                                pe2.profile, subj2, slen2, config.gapopen, config.gapext);
                            hits[clamp_idx].alnscore = result2->score;
                            hits[clamp_idx].qend = static_cast<uint32_t>(result2->end_query);
                            hits[clamp_idx].send = new_ext_start + static_cast<uint32_t>(result2->end_ref);
                            parasail_result_free(result2);
                        }

                        changed = true;
                    }

                    if (changed) {
                        group.indices.erase(
                            std::remove_if(group.indices.begin(), group.indices.end(),
                                [&hit_valid](size_t i) { return !hit_valid[i]; }),
                            group.indices.end());
                        std::sort(group.indices.begin(), group.indices.end(),
                            [&hits](size_t a, size_t b) {
                                return hits[a].sstart < hits[b].sstart;
                            });
                    }
                }
            }
        }

        // Filter survivors into `filtered` and release this batch's
        // batch-local heap.
        for (size_t gi : batch.group_idxs) {
            for (size_t idx : groups[gi].indices) {
                if (!hit_valid[idx] || subject_subseqs[idx].empty()) continue;

                bool keep = true;
                if (config.traceback) {
                    if (config.min_ppositive > 0 && hits[idx].ppositive < config.min_ppositive) keep = false;
                    if (config.min_npositive > 0 && hits[idx].npositive < config.min_npositive) keep = false;
                }
                if (keep) {
                    filtered.push_back(std::move(hits[idx]));
                }
            }
        }

        // Release decoded subseqs and any aligned strings left on dropped
        // hits.  std::string::clear() keeps capacity; shrink_to_fit returns
        // the heap.
        for (size_t gi : batch.group_idxs) {
            for (size_t idx : groups[gi].indices) {
                if (idx < subject_subseqs.size()) {
                    std::string().swap(subject_subseqs[idx]);
                }
                if (idx < hits.size()) {
                    std::string().swap(hits[idx].qseq);
                    std::string().swap(hits[idx].sseq);
                }
            }
        }
    }

    // Free profiles
    for (auto& [key, pe] : profiles) {
        if (pe.profile) parasail_profile_free(pe.profile);
    }

    return filtered;
}

} // namespace ikafssn
