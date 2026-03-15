#include "search/stage3_alignment.hpp"
#include "io/blastdb_reader.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <string>
#include <unordered_map>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/task_arena.h>

#include <parasail.h>

namespace ikafssn {

// Reverse complement a DNA string (returns a new string).
static std::string reverse_complement(const std::string& seq) {
    std::string rc(seq.rbegin(), seq.rend());
    for (auto& c : rc) {
        switch (c) {
            case 'A': c = 'T'; break;
            case 'T': c = 'A'; break;
            case 'C': c = 'G'; break;
            case 'G': c = 'C'; break;
            case 'a': c = 't'; break;
            case 't': c = 'a'; break;
            case 'c': c = 'g'; break;
            case 'g': c = 'c'; break;
            // N and others stay as-is
        }
    }
    return rc;
}

// Walk a parasail CIGAR to compute nident, nmismatch, and build a CIGAR string.
struct CigarStats {
    uint32_t nident = 0;
    uint32_t nmismatch = 0;
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
            case '=': stats.nident += len; break;
            case 'X': stats.nmismatch += len; break;
            case 'I': break; // insertion in query
            case 'D': break; // deletion in query (insertion in ref)
            default: break;
        }
        stats.cigar_str += std::to_string(len);
        stats.cigar_str += op;
    }
    return stats;
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

    // 1. Open BLAST DB volumes
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

    // 2. Build query lookup: query_id -> index in queries[]
    std::unordered_map<std::string, size_t> query_map;
    for (size_t i = 0; i < queries.size(); i++) {
        query_map[queries[i].id] = i;
    }

    // 3. Pre-fetch subject subsequences (volume-parallel)
    // Group hits by volume index (using OutputHit.volume directly)
    std::vector<std::vector<size_t>> hits_by_reader(readers.size());
    std::vector<bool> hit_valid(hits.size(), true);
    for (size_t i = 0; i < hits.size(); i++) {
        if (hits[i].volume < readers.size()) {
            hits_by_reader[hits[i].volume].push_back(i);
        } else {
            logger.warn("Stage 3: hit volume %u out of range (max %zu), skipping",
                        static_cast<unsigned>(hits[i].volume), readers.size());
            hit_valid[i] = false;
        }
    }

    // Sort each volume's hit indices by OID for sequential mmap access
    for (auto& vol_hits : hits_by_reader) {
        std::sort(vol_hits.begin(), vol_hits.end(),
            [&hits](size_t a, size_t b) { return hits[a].oid < hits[b].oid; });
    }

    std::vector<std::string> subject_subseqs(hits.size());
    std::vector<uint32_t> ext_starts(hits.size(), 0);

    int actual_fetch_threads = std::min(config.fetch_threads,
                                         static_cast<int>(readers.size()));
    if (actual_fetch_threads < 1) actual_fetch_threads = 1;

    tbb::task_arena fetch_arena(actual_fetch_threads);
    fetch_arena.execute([&] {
        tbb::parallel_for(size_t(0), readers.size(), [&](size_t ri) {
            for (size_t hit_idx : hits_by_reader[ri]) {
                uint32_t oid = hits[hit_idx].oid;
                uint32_t seq_len = readers[ri].seq_length(oid);

                // Find query length for ratio context
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

                subject_subseqs[hit_idx] = readers[ri].get_subsequence(oid, ext_start, ext_end);
                ext_starts[hit_idx] = ext_start;
                hits[hit_idx].slen = seq_len;
            }
        });
    });

    // 4. Build query profiles (per unique query_id x strand)
    // Key: (query_idx, is_reverse)
    const parasail_matrix_t* matrix = parasail_matrix_lookup("nuc44");

    struct ProfileEntry {
        parasail_profile_t* profile = nullptr;
        std::string seq; // keep alive for profile lifetime
    };
    // We don't share profiles across threads for parasail (profiles are read-only, safe to share)
    std::unordered_map<std::string, ProfileEntry> profiles; // key: "qidx:strand"

    // Collect unique (query_idx, strand) pairs needed
    for (size_t i = 0; i < hits.size(); i++) {
        if (!hit_valid[i] || subject_subseqs[i].empty()) continue;
        auto qit = query_map.find(hits[i].qseqid);
        if (qit == query_map.end()) continue;
        size_t qi = qit->second;
        bool is_rev = (hits[i].sstrand == '-');
        std::string key = std::to_string(qi) + ":" + (is_rev ? "1" : "0");
        if (profiles.find(key) == profiles.end()) {
            ProfileEntry pe;
            if (is_rev) {
                pe.seq = reverse_complement(queries[qi].sequence);
            } else {
                pe.seq = queries[qi].sequence;
            }
            if (config.traceback) {
                pe.profile = parasail_profile_create_sat(
                    pe.seq.c_str(), static_cast<int>(pe.seq.size()), matrix);
            } else {
                pe.profile = parasail_profile_create_sat(
                    pe.seq.c_str(), static_cast<int>(pe.seq.size()), matrix);
            }
            profiles[key] = std::move(pe);
        }
    }

    // 5. Parallel alignment
    // Collect valid hit indices for parallel_for
    std::vector<size_t> valid_indices;
    valid_indices.reserve(hits.size());
    for (size_t i = 0; i < hits.size(); i++) {
        if (hit_valid[i] && !subject_subseqs[i].empty()) {
            auto qit = query_map.find(hits[i].qseqid);
            if (qit != query_map.end()) {
                valid_indices.push_back(i);
            }
        }
    }

    logger.debug("Stage 3: aligning %zu hits (%zu profiles)", valid_indices.size(), profiles.size());

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
            // Traceback alignment
            parasail_result_t* result = parasail_sg_trace_striped_profile_sat(
                pe.profile, subj, slen, config.gapopen, config.gapext);

            hits[idx].alnscore = result->score;

            // Get CIGAR
            parasail_cigar_t* cigar = parasail_result_get_cigar(
                result, pe.seq.c_str(), static_cast<int>(pe.seq.size()),
                subj, slen, matrix);

            // Update coordinates from traceback
            hits[idx].qstart = static_cast<uint32_t>(cigar->beg_query);
            hits[idx].qend = static_cast<uint32_t>(result->end_query);
            hits[idx].sstart = ext_starts[idx] + static_cast<uint32_t>(cigar->beg_ref);
            hits[idx].send = ext_starts[idx] + static_cast<uint32_t>(result->end_ref);

            // Walk CIGAR for stats
            CigarStats cs = walk_cigar(cigar);
            hits[idx].nident = cs.nident;
            hits[idx].mismatch = cs.nmismatch;
            hits[idx].cigar = cs.cigar_str;
            hits[idx].pident = (cs.aln_len > 0) ? 100.0 * cs.nident / cs.aln_len : 0.0;

            // Get traceback strings
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
            // Score-only alignment (no traceback)
            parasail_result_t* result = parasail_sg_striped_profile_sat(
                pe.profile, subj, slen, config.gapopen, config.gapext);

            hits[idx].alnscore = result->score;
            // Update end coordinates from alignment
            hits[idx].qend = static_cast<uint32_t>(result->end_query);
            hits[idx].send = ext_starts[idx] + static_cast<uint32_t>(result->end_ref);
            // q_start and s_start remain from Stage 2 (approximate)

            parasail_result_free(result);
        }
    });

    // 5.5. Overlap resolution for multi-chain hits (context > 0 only)
    // When multiple chains exist for the same (qseqid, sseqid, sstrand),
    // context extension may cause overlapping alignment regions.
    // Clamp the lower-scoring hit's context and re-align.
    {
        bool has_context = context_is_ratio ? (context_ratio > 0) : (context_abs > 0);
        if (has_context) {
            // Group hits by (qseqid, sseqid, sstrand)
            struct HitGroup {
                std::vector<size_t> indices;
            };
            std::unordered_map<std::string, HitGroup> groups;
            for (size_t i = 0; i < hits.size(); i++) {
                if (!hit_valid[i] || subject_subseqs[i].empty()) continue;
                std::string key = hits[i].qseqid + "\t" + hits[i].sseqid + "\t" + hits[i].sstrand;
                groups[key].indices.push_back(i);
            }

            for (auto& [key, group] : groups) {
                if (group.indices.size() < 2) continue;

                // Sort by sstart ascending
                std::sort(group.indices.begin(), group.indices.end(),
                    [&hits](size_t a, size_t b) {
                        return hits[a].sstart < hits[b].sstart;
                    });

                // Iterative overlap resolution
                bool changed = true;
                while (changed) {
                    changed = false;

                    for (size_t gi = 0; gi + 1 < group.indices.size(); gi++) {
                        size_t idx_a = group.indices[gi];
                        size_t idx_b = group.indices[gi + 1];

                        if (!hit_valid[idx_a] || !hit_valid[idx_b]) continue;

                        // Check overlap: hit_a.send >= hit_b.sstart
                        if (hits[idx_a].send < hits[idx_b].sstart) continue;

                        // Determine which has higher chainscore (keep that one intact)
                        size_t keep_idx, clamp_idx;
                        if (hits[idx_a].chainscore >= hits[idx_b].chainscore) {
                            keep_idx = idx_a;
                            clamp_idx = idx_b;
                        } else {
                            keep_idx = idx_b;
                            clamp_idx = idx_a;
                        }

                        // Clamp the lower-scoring hit's overlapping side
                        // If clamp_idx is before keep_idx: clamp its send
                        // If clamp_idx is after keep_idx: clamp its sstart
                        uint32_t new_ext_start, new_ext_end;
                        uint32_t oid = hits[clamp_idx].oid;
                        uint32_t seq_len = hits[clamp_idx].slen;
                        if (seq_len == 0) seq_len = 1; // safety

                        auto qit2 = query_map.find(hits[clamp_idx].qseqid);
                        if (qit2 == query_map.end()) continue;
                        uint32_t query_len2 = static_cast<uint32_t>(queries[qit2->second].sequence.size());
                        uint32_t ctx2 = context_is_ratio
                            ? static_cast<uint32_t>(query_len2 * context_ratio)
                            : context_abs;

                        if (hits[clamp_idx].sstart <= hits[keep_idx].sstart) {
                            // clamp_idx is before keep_idx; clamp its end
                            uint32_t boundary = hits[keep_idx].sstart;
                            // Original stage2 region for clamp hit
                            // We need the stage2 s_end to check if the entire region is consumed
                            if (hits[clamp_idx].sstart >= boundary) {
                                // Entire clamp hit is consumed by keep hit
                                hit_valid[clamp_idx] = false;
                                changed = true;
                                continue;
                            }
                            new_ext_start = ext_starts[clamp_idx];
                            new_ext_end = (boundary > 0) ? boundary - 1 : 0;
                        } else {
                            // clamp_idx is after keep_idx; clamp its start
                            uint32_t boundary = hits[keep_idx].send + 1;
                            if (boundary >= hits[clamp_idx].send) {
                                // Entire clamp hit is consumed
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

                        // Re-fetch subject subsequence
                        uint16_t vol = hits[clamp_idx].volume;
                        if (vol >= readers.size()) {
                            hit_valid[clamp_idx] = false;
                            changed = true;
                            continue;
                        }
                        subject_subseqs[clamp_idx] = readers[vol].get_subsequence(
                            oid, new_ext_start, new_ext_end);
                        ext_starts[clamp_idx] = new_ext_start;

                        // Re-align
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
                            hits[clamp_idx].nident = cs2.nident;
                            hits[clamp_idx].mismatch = cs2.nmismatch;
                            hits[clamp_idx].cigar = cs2.cigar_str;
                            hits[clamp_idx].pident = (cs2.aln_len > 0) ? 100.0 * cs2.nident / cs2.aln_len : 0.0;

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

                    // Re-sort by sstart after changes
                    if (changed) {
                        // Remove invalidated indices
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
    }

    // 6. Filter by min_pident / min_nident (only meaningful with traceback)
    std::vector<OutputHit> filtered;
    filtered.reserve(hits.size());
    for (size_t i = 0; i < hits.size(); i++) {
        if (!hit_valid[i] || subject_subseqs[i].empty()) continue;

        // Check query exists
        auto qit = query_map.find(hits[i].qseqid);
        if (qit == query_map.end()) continue;

        if (config.traceback) {
            if (config.min_pident > 0 && hits[i].pident < config.min_pident) continue;
            if (config.min_nident > 0 && hits[i].nident < config.min_nident) continue;
        }
        filtered.push_back(std::move(hits[i]));
    }

    // 7. Free profiles
    for (auto& [key, pe] : profiles) {
        if (pe.profile) parasail_profile_free(pe.profile);
    }

    return filtered;
}

} // namespace ikafssn
