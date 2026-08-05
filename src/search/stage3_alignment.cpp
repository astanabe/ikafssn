#include "search/stage3_alignment.hpp"
#include "core/spaced_seed.hpp"
#include "io/blastdb_reader.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#include <parasail.h>

namespace ikafssn {

// Walk a parasail CIGAR and report the aligned region only.
//
// A semi-global CIGAR materialises the free end gaps: the traceback reaches
// the origin, so `cigar->beg_query` / `beg_ref` are always 0, and the tail gap
// of whichever sequence has not reached its end is appended.  The run
// therefore starts and ends with gap-only ops describing the context window
// rather than the match, and everything outside the first and last '=' / 'X'
// is stripped.  The ends need no such trim: `result->end_query` / `end_ref`
// are the traceback's starting cell and already sit on the last aligned base.
struct CigarStats {
    uint32_t npositive = 0;
    uint32_t nnegative = 0;
    uint32_t aln_len = 0;      // columns in the aligned region
    uint32_t lead_query = 0;   // query bases before the first aligned column
    uint32_t lead_subject = 0; // subject bases before the first aligned column
    uint32_t lead_cols = 0;    // columns before the first aligned column
    std::string cigar_str;     // aligned region only
};

// 'I' consumes the query only and 'D' the subject only, following SAM.
static CigarStats walk_cigar(const parasail_cigar_t* cigar) {
    int first = -1, last = -1;
    for (int i = 0; i < cigar->len; i++) {
        char op = parasail_cigar_decode_op(cigar->seq[i]);
        if (op == '=' || op == 'X') {
            if (first < 0) first = i;
            last = i;
        }
    }
    // No aligned column: nothing to anchor the trim to, so keep the run as is.
    if (first < 0) { first = 0; last = cigar->len - 1; }

    CigarStats stats;
    for (int i = 0; i < first; i++) {
        char op = parasail_cigar_decode_op(cigar->seq[i]);
        uint32_t len = parasail_cigar_decode_len(cigar->seq[i]);
        stats.lead_cols += len;
        if (op == 'I') stats.lead_query += len;
        else if (op == 'D') stats.lead_subject += len;
    }
    for (int i = first; i <= last; i++) {
        char op = parasail_cigar_decode_op(cigar->seq[i]);
        uint32_t len = parasail_cigar_decode_len(cigar->seq[i]);
        stats.aln_len += len;
        if (op == '=') stats.npositive += len;
        else if (op == 'X') stats.nnegative += len;
        stats.cigar_str += std::to_string(len);
        stats.cigar_str += op;
    }
    return stats;
}

// Slice the aligned region out of one parasail traceback row, which holds one
// character per CIGAR column.
static std::string aligned_slice(const char* row, const CigarStats& cs) {
    return std::string(row + cs.lead_cols, cs.aln_len);
}

// ---------------------------------------------------------------------------
// Stage 3 batching
//
// Hits sharing the same (qseqid, sseqid, sstrand) form an atomic group: the
// overlap-resolution loop assumes all such hits are visible at once, so a
// group cannot be split across batches.  Groups are bin-packed greedily
// against `posting_budget` (the residual heap budget plumbed from
// `-memory_limit`).  A group whose own cost exceeds the budget is processed
// as its own solo batch.
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
    uint64_t subseq = static_cast<uint64_t>(h.send - h.sstart) + 2ull * ctx;
    uint64_t tb = cfg.traceback
        ? 2ull * std::max<uint64_t>(h.qlen, subseq)
        : 0ull;
    return subseq + tb + 256; // const overhead: OutputHit + std::string headers + pad
}

// One atomic group, stored as a slice of the flat `group_hits` buffer that
// holds the hit indices of every group back to back.
struct Stage3Group {
    uint32_t begin = 0;      // offset into group_hits
    uint32_t len = 0;        // live hits, shrunk by the overlap resolution
    uint32_t owned_len = 0;  // hits the group owns, dropped ones included
    uint64_t cost = 0;
};

struct Stage3Batch {
    std::vector<size_t> group_idxs; // indices into the ordered groups vector
    uint64_t cost = 0;
    bool oversize = false;
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
            solo.oversize = true;
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

    return run_stage3(hits, queries, readers, config,
                      context_is_ratio, context_ratio, context_abs, logger);
}

std::vector<OutputHit> run_stage3(
    std::vector<OutputHit>& hits,
    const std::vector<FastaRecord>& queries,
    const std::vector<BlastDbReader>& readers,
    const Stage3Config& config,
    bool context_is_ratio,
    double context_ratio,
    uint32_t context_abs,
    const Logger& logger)
{
    if (hits.empty()) return {};

    if (readers.empty()) {
        logger.error("Stage 3: no BLAST DB volumes are open");
        return {};
    }

    // Reject the hits that cannot be aligned: a volume outside `readers` or a
    // query_idx outside `queries`.  Every hit valid from here on indexes
    // `queries` directly through its query_idx.
    std::vector<bool> hit_valid(hits.size(), true);
    for (size_t i = 0; i < hits.size(); i++) {
        if (hits[i].volume >= readers.size()) {
            logger.warn("Stage 3: hit volume %u out of range (max %zu), skipping",
                        static_cast<unsigned>(hits[i].volume), readers.size());
            hit_valid[i] = false;
            continue;
        }
        if (hits[i].query_idx >= queries.size()) {
            logger.warn("Stage 3: hit query index %u out of range (max %zu), skipping",
                        static_cast<unsigned>(hits[i].query_idx), queries.size());
            hit_valid[i] = false;
        }
    }

    // Build profiles up front (per unique query_idx x strand).  Profiles are
    // small and shared across all batches.
    const parasail_matrix_t* matrix = parasail_matrix_lookup(config.score_matrix.c_str());
    if (!matrix) {
        logger.error("Stage 3: unknown score matrix '%s'", config.score_matrix.c_str());
        return {};
    }

    struct ProfileEntry {
        parasail_profile_t* profile = nullptr;
        std::string seq; // keep alive for profile lifetime
        bool built = false;
    };
    // Slot `query_idx * 2 + is_rev`; only the (query, strand) pairs the hits
    // actually use are built.
    std::vector<ProfileEntry> profiles(queries.size() * 2);
    auto profile_slot = [&](size_t hit_idx) {
        return static_cast<size_t>(hits[hit_idx].query_idx) * 2
             + (hits[hit_idx].sstrand == '-' ? 1u : 0u);
    };

    for (size_t i = 0; i < hits.size(); i++) {
        if (!hit_valid[i]) continue;
        ProfileEntry& pe = profiles[profile_slot(i)];
        if (pe.built) continue;
        const size_t qi = hits[i].query_idx;
        pe.seq = (hits[i].sstrand == '-')
            ? reverse_complement_string(queries[qi].sequence)
            : queries[qi].sequence;
        pe.profile = parasail_profile_create_sat(
            pe.seq.c_str(), static_cast<int>(pe.seq.size()), matrix);
        pe.built = true;
    }

    // Build atomic groups keyed by (query index, sseqid, sstrand).  The
    // overlap-resolution loop assumes each group is visible whole, so groups
    // are the unit of batching; every group's hit indices live in one flat
    // buffer, each group owning the slice [begin, begin + len).
    //
    // Sorting by the group key gives contiguous runs, and relabelling those
    // runs while walking the hits in ascending index order puts the groups in
    // first-appearance order with each slice in arrival order.
    std::vector<uint32_t> valid_hits;
    valid_hits.reserve(hits.size());
    for (size_t i = 0; i < hits.size(); i++)
        if (hit_valid[i]) valid_hits.push_back(static_cast<uint32_t>(i));

    std::vector<uint32_t> run_of(hits.size(), 0);
    uint32_t nruns = 0;
    {
        auto same_key = [&](uint32_t a, uint32_t b) {
            return hits[a].query_idx == hits[b].query_idx
                && hits[a].sstrand == hits[b].sstrand
                && hits[a].sseqid  == hits[b].sseqid;
        };
        std::vector<uint32_t> by_key = valid_hits;
        std::sort(by_key.begin(), by_key.end(), [&](uint32_t a, uint32_t b) {
            if (hits[a].query_idx != hits[b].query_idx)
                return hits[a].query_idx < hits[b].query_idx;
            if (hits[a].sstrand != hits[b].sstrand)
                return hits[a].sstrand < hits[b].sstrand;
            return hits[a].sseqid < hits[b].sseqid;
        });
        for (size_t i = 0; i < by_key.size(); i++) {
            if (i != 0 && !same_key(by_key[i - 1], by_key[i])) ++nruns;
            run_of[by_key[i]] = nruns;
        }
        if (!by_key.empty()) ++nruns;
    }

    std::vector<uint32_t> group_of_run(nruns, UINT32_MAX);
    std::vector<Stage3Group> groups;
    for (uint32_t h : valid_hits) {
        uint32_t& g = group_of_run[run_of[h]];
        if (g == UINT32_MAX) {
            g = static_cast<uint32_t>(groups.size());
            groups.push_back(Stage3Group{});
        }
        ++groups[g].owned_len;
    }
    uint32_t acc = 0;
    for (auto& g : groups) {
        g.begin = acc;
        acc += g.owned_len;
    }
    std::vector<uint32_t> group_hits(valid_hits.size());
    for (uint32_t h : valid_hits) {
        Stage3Group& g = groups[group_of_run[run_of[h]]];
        group_hits[g.begin + g.len++] = h;
    }

    for (auto& g : groups) {
        // Sort by sstart for the overlap-resolution loop's invariant.
        auto slice = group_hits.begin() + g.begin;
        std::sort(slice, slice + g.len,
            [&hits](uint32_t a, uint32_t b) { return hits[a].sstart < hits[b].sstart; });
        for (uint32_t t = 0; t < g.len; t++) {
            g.cost += compute_hit_cost(hits[group_hits[g.begin + t]], config,
                                       context_is_ratio, context_ratio, context_abs);
        }
    }

    // Plan the batches.
    auto batches = plan_stage3_batches(groups, config.posting_budget);

    {
        size_t oversize_count = 0;
        for (const auto& b : batches) if (b.oversize) ++oversize_count;
        logger.info("Stage 3 batch plan: %zu batch(es) over %zu group(s) "
                    "(budget=%llu, oversize=%zu)",
                    batches.size(), groups.size(),
                    static_cast<unsigned long long>(config.posting_budget),
                    oversize_count);
    }

    // Per-hit scratch storage.  Sized to hits.size() so we can index by the
    // hit index, but each batch only fills entries for its own hits and
    // clears them when the batch ends.
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
            const auto& g = groups[gi];
            for (uint32_t t = 0; t < g.len; t++) {
                size_t hidx = group_hits[g.begin + t];
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
                     batch.oversize ? " [oversize]" : "");
        if (batch.oversize) {
            logger.warn("Stage 3 batch %zu/%zu: single group exceeds posting_budget "
                        "(%llu > %llu); processing solo (over-budget)",
                        bi + 1, batches.size(),
                        static_cast<unsigned long long>(batch.cost),
                        static_cast<unsigned long long>(config.posting_budget));
        }

        // Fetch subseqs in (volume, OID) order via a single hit-parallel
        // parallel_for (default arena = -nthread), preserving sequential mmap
        // locality within each task.  `get_subsequence` is lock-free (no
        // CSeqDBAtlas lock; mmap arithmetic + ncbi2na/ambig decode only), so
        // concurrent calls are safe.
        std::vector<size_t> ordered_hits;
        ordered_hits.reserve(batch_hit_count);
        for (size_t ri = 0; ri < readers.size(); ri++) {
            for (size_t h : hits_by_reader[ri]) ordered_hits.push_back(h);
        }
        // The hint is process-wide, so one call covers every open volume.
        readers.front().set_mmap_strategy(BlastDbReader::MMapStrategy::kNormal);
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, ordered_hits.size(), 16),
            [&](const tbb::blocked_range<size_t>& r) {
                for (size_t i = r.begin(); i < r.end(); i++) {
                    size_t hit_idx = ordered_hits[i];
                    uint16_t vol = hits[hit_idx].volume;
                    uint32_t oid = hits[hit_idx].oid;

                    uint32_t query_len = static_cast<uint32_t>(
                        queries[hits[hit_idx].query_idx].sequence.size());

                    uint32_t ctx = context_is_ratio
                        ? static_cast<uint32_t>(query_len * context_ratio)
                        : context_abs;

                    ContextSubseq cs = extract_context_subseq(
                        readers[vol], oid, hits[hit_idx].sstart, hits[hit_idx].send, ctx);
                    subject_subseqs[hit_idx] = std::move(cs.seq);
                    ext_starts[hit_idx] = cs.ext_start;
                    hits[hit_idx].slen = cs.seq_len;
                }
            });
        // Subseqs are now in `subject_subseqs` (heap-owned strings); the
        // alignment step does not touch the mmap further this batch, so
        // release the cached pages before the next batch's fetch.
        readers.front().set_mmap_strategy(BlastDbReader::MMapStrategy::kDontNeed);

        // Collect valid hit indices for parallel alignment.
        std::vector<size_t> valid_indices;
        valid_indices.reserve(batch_hit_count);
        for (size_t gi : batch.group_idxs) {
            const auto& g = groups[gi];
            for (uint32_t t = 0; t < g.len; t++) {
                size_t hidx = group_hits[g.begin + t];
                if (hit_valid[hidx] && !subject_subseqs[hidx].empty()) {
                    valid_indices.push_back(hidx);
                }
            }
        }

        // Align this batch in parallel.
        tbb::parallel_for(size_t(0), valid_indices.size(), [&](size_t vi) {
            size_t idx = valid_indices[vi];
            const ProfileEntry& pe = profiles[profile_slot(idx)];

            const char* subj = subject_subseqs[idx].c_str();
            int slen = static_cast<int>(subject_subseqs[idx].size());

            if (config.traceback) {
                parasail_result_t* result = parasail_sg_trace_striped_profile_sat(
                    pe.profile, subj, slen, config.gapopen, config.gapext);

                hits[idx].alnscore = result->score;

                parasail_cigar_t* cigar = parasail_result_get_cigar(
                    result, pe.seq.c_str(), static_cast<int>(pe.seq.size()),
                    subj, slen, matrix);
                if (!cigar) {
                    parasail_result_free(result);
                    return;
                }

                CigarStats cs = walk_cigar(cigar);
                hits[idx].qstart = static_cast<uint32_t>(cigar->beg_query) + cs.lead_query;
                hits[idx].qend = static_cast<uint32_t>(result->end_query) + 1;
                hits[idx].sstart = ext_starts[idx] +
                                   static_cast<uint32_t>(cigar->beg_ref) + cs.lead_subject;
                hits[idx].send = ext_starts[idx] +
                                 static_cast<uint32_t>(result->end_ref) + 1;

                hits[idx].npositive = cs.npositive;
                hits[idx].nnegative = cs.nnegative;
                hits[idx].cigar = cs.cigar_str;
                hits[idx].ppositive = (cs.aln_len > 0) ? 100.0 * cs.npositive / cs.aln_len : 0.0;

                parasail_traceback_t* tb = parasail_result_get_traceback(
                    result, pe.seq.c_str(), static_cast<int>(pe.seq.size()),
                    subj, slen, matrix, '|', '*', ' ');
                if (tb) {
                    hits[idx].qseq = aligned_slice(tb->query, cs);
                    hits[idx].sseq = aligned_slice(tb->ref, cs);
                    parasail_traceback_free(tb);
                }

                parasail_cigar_free(cigar);
                parasail_result_free(result);
            } else {
                parasail_result_t* result = parasail_sg_striped_profile_sat(
                    pe.profile, subj, slen, config.gapopen, config.gapext);

                hits[idx].alnscore = result->score;
                hits[idx].qend = static_cast<uint32_t>(result->end_query) + 1;
                hits[idx].send = ext_starts[idx] +
                                 static_cast<uint32_t>(result->end_ref) + 1;

                parasail_result_free(result);
            }
        });

        // Overlap resolution for multi-chain hits in this batch's groups
        // (context > 0 only).  Groups in different batches share no state,
        // so the loop is naturally batch-local.
        if (has_context) {
            for (size_t gi : batch.group_idxs) {
                auto& group = groups[gi];
                if (group.len < 2) continue;

                bool changed = true;
                while (changed) {
                    changed = false;

                    for (uint32_t pi = 0; pi + 1 < group.len; pi++) {
                        size_t idx_a = group_hits[group.begin + pi];
                        size_t idx_b = group_hits[group.begin + pi + 1];

                        if (!hit_valid[idx_a] || !hit_valid[idx_b]) continue;
                        if (hits[idx_a].send <= hits[idx_b].sstart) continue;

                        size_t keep_idx, clamp_idx;
                        if (hits[idx_a].chainscore >= hits[idx_b].chainscore) {
                            keep_idx = idx_a;
                            clamp_idx = idx_b;
                        } else {
                            keep_idx = idx_b;
                            clamp_idx = idx_a;
                        }

                        // Half-open, like the hit coordinates it derives from.
                        uint32_t new_ext_start, new_ext_end;
                        uint32_t oid = hits[clamp_idx].oid;
                        uint32_t seq_len = hits[clamp_idx].slen;
                        if (seq_len == 0) seq_len = 1;

                        uint32_t query_len2 = static_cast<uint32_t>(
                            queries[hits[clamp_idx].query_idx].sequence.size());
                        uint32_t ctx2 = context_is_ratio
                            ? static_cast<uint32_t>(query_len2 * context_ratio)
                            : context_abs;

                        if (hits[clamp_idx].sstart <= hits[keep_idx].sstart) {
                            // Clamp the left hit to end where the kept one starts.
                            uint32_t boundary = hits[keep_idx].sstart;
                            if (hits[clamp_idx].sstart >= boundary) {
                                hit_valid[clamp_idx] = false;
                                changed = true;
                                continue;
                            }
                            new_ext_start = ext_starts[clamp_idx];
                            new_ext_end = boundary;
                        } else {
                            // Clamp the right hit to start where the kept one ends.
                            uint32_t boundary = hits[keep_idx].send;
                            if (boundary >= hits[clamp_idx].send) {
                                hit_valid[clamp_idx] = false;
                                changed = true;
                                continue;
                            }
                            new_ext_start = boundary;
                            new_ext_end = std::min<uint64_t>(
                                static_cast<uint64_t>(hits[clamp_idx].send) + ctx2, seq_len);
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
                        // get_subsequence takes an inclusive end.
                        subject_subseqs[clamp_idx] = readers[vol].get_subsequence(
                            oid, new_ext_start, new_ext_end - 1);
                        ext_starts[clamp_idx] = new_ext_start;

                        const ProfileEntry& pe2 = profiles[profile_slot(clamp_idx)];

                        const char* subj2 = subject_subseqs[clamp_idx].c_str();
                        int slen2 = static_cast<int>(subject_subseqs[clamp_idx].size());

                        if (config.traceback) {
                            parasail_result_t* result2 = parasail_sg_trace_striped_profile_sat(
                                pe2.profile, subj2, slen2, config.gapopen, config.gapext);
                            hits[clamp_idx].alnscore = result2->score;

                            parasail_cigar_t* cigar2 = parasail_result_get_cigar(
                                result2, pe2.seq.c_str(), static_cast<int>(pe2.seq.size()),
                                subj2, slen2, matrix);
                            if (!cigar2) {
                                hit_valid[clamp_idx] = false;
                                parasail_result_free(result2);
                                changed = true;
                                continue;
                            }

                            CigarStats cs2 = walk_cigar(cigar2);
                            hits[clamp_idx].qstart =
                                static_cast<uint32_t>(cigar2->beg_query) + cs2.lead_query;
                            hits[clamp_idx].qend =
                                static_cast<uint32_t>(result2->end_query) + 1;
                            hits[clamp_idx].sstart = new_ext_start +
                                static_cast<uint32_t>(cigar2->beg_ref) + cs2.lead_subject;
                            hits[clamp_idx].send = new_ext_start +
                                static_cast<uint32_t>(result2->end_ref) + 1;

                            hits[clamp_idx].npositive = cs2.npositive;
                            hits[clamp_idx].nnegative = cs2.nnegative;
                            hits[clamp_idx].cigar = cs2.cigar_str;
                            hits[clamp_idx].ppositive = (cs2.aln_len > 0) ? 100.0 * cs2.npositive / cs2.aln_len : 0.0;

                            parasail_traceback_t* tb2 = parasail_result_get_traceback(
                                result2, pe2.seq.c_str(), static_cast<int>(pe2.seq.size()),
                                subj2, slen2, matrix, '|', '*', ' ');
                            if (tb2) {
                                hits[clamp_idx].qseq = aligned_slice(tb2->query, cs2);
                                hits[clamp_idx].sseq = aligned_slice(tb2->ref, cs2);
                                parasail_traceback_free(tb2);
                            }

                            parasail_cigar_free(cigar2);
                            parasail_result_free(result2);
                        } else {
                            parasail_result_t* result2 = parasail_sg_striped_profile_sat(
                                pe2.profile, subj2, slen2, config.gapopen, config.gapext);
                            hits[clamp_idx].alnscore = result2->score;
                            hits[clamp_idx].qend =
                                static_cast<uint32_t>(result2->end_query) + 1;
                            hits[clamp_idx].send = new_ext_start +
                                static_cast<uint32_t>(result2->end_ref) + 1;
                            parasail_result_free(result2);
                        }

                        changed = true;
                    }

                    if (changed) {
                        // Partition rather than erase: the dropped hits stay
                        // in the slice's tail so the batch's release loop
                        // still sees them.
                        auto slice = group_hits.begin() + group.begin;
                        auto live_end = std::stable_partition(
                            slice, slice + group.len,
                            [&hit_valid](uint32_t i) { return hit_valid[i]; });
                        group.len = static_cast<uint32_t>(live_end - slice);
                        std::sort(slice, slice + group.len,
                            [&hits](uint32_t a, uint32_t b) {
                                return hits[a].sstart < hits[b].sstart;
                            });
                    }
                }
            }
        }

        // Filter the survivors into `filtered`.
        for (size_t gi : batch.group_idxs) {
            const auto& g = groups[gi];
            for (uint32_t t = 0; t < g.len; t++) {
                size_t idx = group_hits[g.begin + t];
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

        // Release this batch's decoded subseqs and aligned strings.  Walks
        // the whole slice so the hits the overlap resolution dropped hand
        // their heap back too; swapping in an empty string returns the
        // buffer, which clear() would keep.
        for (size_t gi : batch.group_idxs) {
            const auto& g = groups[gi];
            for (uint32_t t = 0; t < g.owned_len; t++) {
                size_t idx = group_hits[g.begin + t];
                std::string().swap(subject_subseqs[idx]);
                std::string().swap(hits[idx].qseq);
                std::string().swap(hits[idx].sseq);
            }
        }
    }

    // Free the profiles.
    for (auto& pe : profiles) {
        if (pe.profile) parasail_profile_free(pe.profile);
    }

    return filtered;
}

} // namespace ikafssn
