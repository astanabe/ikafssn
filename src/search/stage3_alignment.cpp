#include "search/stage3_alignment.hpp"
#include "io/blastdb_reader.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
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

    // 1. Open BLAST DB volumes and build accession -> (reader_idx, oid) map
    auto vol_paths = BlastDbReader::find_volume_paths(db_path);
    if (vol_paths.empty()) {
        logger.error("Stage 3: no BLAST DB volumes found at '%s'", db_path.c_str());
        return {};
    }

    std::vector<BlastDbReader> readers(vol_paths.size());
    std::unordered_map<std::string, std::pair<size_t, uint32_t>> acc_map;

    for (size_t vi = 0; vi < vol_paths.size(); vi++) {
        if (!readers[vi].open(vol_paths[vi])) {
            logger.error("Stage 3: cannot open volume '%s'", vol_paths[vi].c_str());
            return {};
        }
        uint32_t nseqs = readers[vi].num_sequences();
        for (uint32_t oid = 0; oid < nseqs; oid++) {
            std::string acc = readers[vi].get_accession(oid);
            if (!acc.empty()) {
                acc_map[acc] = {vi, oid};
            }
        }
    }

    // 2. Build query lookup: query_id -> index in queries[]
    std::unordered_map<std::string, size_t> query_map;
    for (size_t i = 0; i < queries.size(); i++) {
        query_map[queries[i].id] = i;
    }

    // 3. Pre-fetch subject subsequences (volume-parallel)
    // Group hits by reader_idx
    std::vector<std::vector<size_t>> hits_by_reader(readers.size());
    std::vector<bool> hit_valid(hits.size(), true);
    for (size_t i = 0; i < hits.size(); i++) {
        auto it = acc_map.find(hits[i].accession);
        if (it != acc_map.end()) {
            hits_by_reader[it->second.first].push_back(i);
        } else {
            logger.warn("Stage 3: accession '%s' not found in BLAST DB, skipping",
                        hits[i].accession.c_str());
            hit_valid[i] = false;
        }
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
                auto it = acc_map.find(hits[hit_idx].accession);
                uint32_t oid = it->second.second;
                uint32_t seq_len = readers[ri].seq_length(oid);

                // Find query length for ratio context
                uint32_t query_len = 0;
                auto qit = query_map.find(hits[hit_idx].query_id);
                if (qit != query_map.end()) {
                    query_len = static_cast<uint32_t>(queries[qit->second].sequence.size());
                }

                uint32_t ctx = context_is_ratio
                    ? static_cast<uint32_t>(query_len * context_ratio)
                    : context_abs;

                uint32_t ext_start = (hits[hit_idx].s_start >= ctx)
                    ? hits[hit_idx].s_start - ctx : 0;
                uint32_t ext_end = std::min(hits[hit_idx].s_end + ctx, seq_len - 1);

                std::string full_seq = readers[ri].get_sequence(oid);
                if (!full_seq.empty() && ext_start <= ext_end && ext_end < full_seq.size()) {
                    subject_subseqs[hit_idx] = full_seq.substr(ext_start, ext_end - ext_start + 1);
                } else if (!full_seq.empty()) {
                    // Clamp
                    if (ext_end >= full_seq.size()) ext_end = static_cast<uint32_t>(full_seq.size()) - 1;
                    if (ext_start <= ext_end) {
                        subject_subseqs[hit_idx] = full_seq.substr(ext_start, ext_end - ext_start + 1);
                    }
                }
                ext_starts[hit_idx] = ext_start;
                hits[hit_idx].s_length = seq_len;
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
        auto qit = query_map.find(hits[i].query_id);
        if (qit == query_map.end()) continue;
        size_t qi = qit->second;
        bool is_rev = (hits[i].strand == '-');
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
            auto qit = query_map.find(hits[i].query_id);
            if (qit != query_map.end()) {
                valid_indices.push_back(i);
            }
        }
    }

    logger.debug("Stage 3: aligning %zu hits (%zu profiles)", valid_indices.size(), profiles.size());

    tbb::parallel_for(size_t(0), valid_indices.size(), [&](size_t vi) {
        size_t idx = valid_indices[vi];
        auto qit = query_map.find(hits[idx].query_id);
        size_t qi = qit->second;
        bool is_rev = (hits[idx].strand == '-');
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
            hits[idx].q_start = static_cast<uint32_t>(cigar->beg_query);
            hits[idx].q_end = static_cast<uint32_t>(result->end_query);
            hits[idx].s_start = ext_starts[idx] + static_cast<uint32_t>(cigar->beg_ref);
            hits[idx].s_end = ext_starts[idx] + static_cast<uint32_t>(result->end_ref);

            // Walk CIGAR for stats
            CigarStats cs = walk_cigar(cigar);
            hits[idx].nident = cs.nident;
            hits[idx].nmismatch = cs.nmismatch;
            hits[idx].cigar = cs.cigar_str;
            hits[idx].pident = (cs.aln_len > 0) ? 100.0 * cs.nident / cs.aln_len : 0.0;

            // Get traceback strings
            parasail_traceback_t* tb = parasail_result_get_traceback(
                result, pe.seq.c_str(), static_cast<int>(pe.seq.size()),
                subj, slen, matrix, '|', '*', ' ');
            if (tb) {
                hits[idx].q_seq = tb->query;
                hits[idx].s_seq = tb->ref;
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
            hits[idx].q_end = static_cast<uint32_t>(result->end_query);
            hits[idx].s_end = ext_starts[idx] + static_cast<uint32_t>(result->end_ref);
            // q_start and s_start remain from Stage 2 (approximate)

            parasail_result_free(result);
        }
    });

    // 6. Filter by min_pident / min_nident (only meaningful with traceback)
    std::vector<OutputHit> filtered;
    filtered.reserve(hits.size());
    for (size_t i = 0; i < hits.size(); i++) {
        if (!hit_valid[i] || subject_subseqs[i].empty()) continue;

        // Check query exists
        auto qit = query_map.find(hits[i].query_id);
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
