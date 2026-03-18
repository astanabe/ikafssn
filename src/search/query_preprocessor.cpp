#include "search/query_preprocessor.hpp"
#include "search/volume_searcher.hpp"
#include "search/stage1_filter.hpp"
#include "index/kix_reader.hpp"
#include "index/khx_reader.hpp"
#include "core/kmer_encoding.hpp"
#include "core/spaced_seed.hpp"
#include "core/config.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <unordered_set>

namespace ikafssn {

template <typename KmerInt>
static std::vector<std::pair<uint32_t, KmerInt>>
extract_kmers(const std::string& seq, int k, bool* has_multi_degen = nullptr,
              int max_expansion = 16) {
    std::vector<std::pair<uint32_t, KmerInt>> kmers;
    KmerScanner<KmerInt> scanner(k);
    scanner.scan_ambig(seq.data(), seq.size(),
        [&](uint32_t pos, KmerInt kmer) {
            kmers.emplace_back(pos, kmer);
        },
        [&](uint32_t pos, KmerInt base_kmer, const AmbigInfo* infos, int count) {
            expand_ambig_kmer_multi<KmerInt>(base_kmer, infos, count,
                [&](KmerInt expanded) {
                    kmers.emplace_back(pos, expanded);
                });
        },
        has_multi_degen,
        max_expansion);
    return kmers;
}

template <typename KmerInt>
static std::vector<std::pair<uint32_t, KmerInt>>
extract_kmers_spaced(const std::string& seq, int k,
                     const std::vector<uint32_t>& masks, int t,
                     const std::vector<KmerInt>& mask_tags = {},
                     bool* has_multi_degen = nullptr,
                     int max_expansion = 16) {
    std::vector<std::pair<uint32_t, KmerInt>> kmers;
    KmerScanner<KmerInt> scanner(k);
    scanner.scan_spaced_ambig(seq.data(), seq.size(), masks, t,
        [&](uint32_t pos, KmerInt kmer) {
            kmers.emplace_back(pos, kmer);
        },
        [&](uint32_t pos, KmerInt base_kmer, const AmbigInfo* infos, int count) {
            expand_ambig_kmer_multi<KmerInt>(base_kmer, infos, count,
                [&](KmerInt expanded) {
                    kmers.emplace_back(pos, expanded);
                });
        },
        has_multi_degen,
        max_expansion,
        mask_tags);
    return kmers;
}

// Compute global max_freq: aggregate across all volumes if auto mode.
static uint32_t compute_global_max_freq(
    uint32_t config_max_freq,
    const std::vector<const KixReader*>& all_kix) {
    if (config_max_freq > 0) return config_max_freq;

    // Auto mode: aggregate total_postings and table_size across all volumes
    uint64_t total_postings = 0;
    uint64_t tbl_size = 0;
    for (const auto* kix : all_kix) {
        total_postings += kix->total_postings();
        tbl_size = kix->table_size(); // same for all volumes
    }
    return compute_effective_max_freq(0, total_postings, tbl_size);
}

template <typename KmerInt>
QueryKmerData<KmerInt> preprocess_query(
    const std::string& query_seq, int k,
    const std::vector<const KixReader*>& all_kix,
    const KhxReader* khx,
    const SearchConfig& config,
    uint8_t t,
    const std::vector<uint32_t>& masks,
    const std::vector<KmerInt>& mask_tags) {

    QueryKmerData<KmerInt> result;

    // 1. Extract forward k-mers
    std::vector<std::pair<uint32_t, KmerInt>> fwd_kmers;
    if (t > 0 && !masks.empty()) {
        fwd_kmers = extract_kmers_spaced<KmerInt>(query_seq, k, masks,
                        static_cast<int>(t), mask_tags, &result.has_multi_degen,
                        static_cast<int>(config.max_degen_expand));
    } else {
        fwd_kmers = extract_kmers<KmerInt>(query_seq, k, &result.has_multi_degen,
                        static_cast<int>(config.max_degen_expand));
    }
    if (fwd_kmers.empty()) return result;

    // 2. Build reverse complement k-mers
    std::vector<std::pair<uint32_t, KmerInt>> rc_kmers;
    if (t > 0 && !masks.empty()) {
        // Spaced seed: scan RC string with same templates, remap positions
        std::string rc_seq = reverse_complement_string(query_seq);
        auto rc_raw = extract_kmers_spaced<KmerInt>(rc_seq, k, masks,
                          static_cast<int>(t), mask_tags, &result.has_multi_degen,
                          static_cast<int>(config.max_degen_expand));
        int span = static_cast<int>(t);
        rc_kmers.reserve(rc_raw.size());
        for (const auto& [pos, kmer] : rc_raw) {
            // Remap position: rc position p corresponds to fwd position (len - p - span)
            uint32_t fwd_pos = static_cast<uint32_t>(query_seq.size()) - pos - static_cast<uint32_t>(span);
            rc_kmers.emplace_back(fwd_pos, kmer);
        }
    } else {
        // Contiguous: existing kmer_revcomp() approach
        rc_kmers.reserve(fwd_kmers.size());
        for (const auto& [pos, kmer] : fwd_kmers) {
            rc_kmers.emplace_back(pos, kmer_revcomp(kmer, k));
        }
    }

    // 3. Determine global max_freq
    uint32_t global_max_freq = compute_global_max_freq(
        config.stage1.max_freq, all_kix);

    // 4-5. Build high-freq set and filter (skip entirely when disabled)
    std::unordered_set<uint32_t> highfreq_set;

    // Helper: convert pair vector to SoA
    auto to_soa = [](const std::vector<std::pair<uint32_t, KmerInt>>& pairs,
                     std::vector<uint32_t>& positions,
                     std::vector<KmerInt>& kmer_values) {
        positions.reserve(pairs.size());
        kmer_values.reserve(pairs.size());
        for (const auto& [pos, kmer] : pairs) {
            positions.push_back(pos);
            kmer_values.push_back(kmer);
        }
    };

    if (global_max_freq == Stage1Config::MAX_FREQ_DISABLED) {
        // High-freq filtering disabled: use all k-mers as-is
        to_soa(fwd_kmers, result.fwd_positions, result.fwd_kmer_values);
        to_soa(rc_kmers, result.rc_positions, result.rc_kmer_values);
    } else {
        // Collect all distinct k-mer values from both strands
        std::unordered_set<uint32_t> all_query_kmer_values;
        for (const auto& [pos, kmer] : fwd_kmers) {
            all_query_kmer_values.insert(static_cast<uint32_t>(kmer));
        }
        for (const auto& [pos, kmer] : rc_kmers) {
            all_query_kmer_values.insert(static_cast<uint32_t>(kmer));
        }

        for (uint32_t kmer_idx : all_query_kmer_values) {
            // Strip mask tag bit for KHX lookup (KHX uses base 4^k table)
            uint32_t khx_idx = kmer_idx & (ikafssn::table_size(k) - 1);
            // Check .khx exclusion
            if (khx != nullptr && khx->is_excluded(khx_idx)) {
                highfreq_set.insert(kmer_idx);
                continue;
            }
            // Sum counts across all volumes
            uint64_t total_count = 0;
            for (const auto* kix : all_kix) {
                total_count += kix->count_postings(kmer_idx);
            }
            if (total_count > global_max_freq) {
                highfreq_set.insert(kmer_idx);
            }
        }

        // Filter out high-freq k-mers from both strand vectors (SoA)
        result.fwd_positions.reserve(fwd_kmers.size());
        result.fwd_kmer_values.reserve(fwd_kmers.size());
        for (const auto& [pos, kmer] : fwd_kmers) {
            if (highfreq_set.count(static_cast<uint32_t>(kmer)) == 0) {
                result.fwd_positions.push_back(pos);
                result.fwd_kmer_values.push_back(kmer);
            }
        }

        result.rc_positions.reserve(rc_kmers.size());
        result.rc_kmer_values.reserve(rc_kmers.size());
        for (const auto& [pos, kmer] : rc_kmers) {
            if (highfreq_set.count(static_cast<uint32_t>(kmer)) == 0) {
                result.rc_positions.push_back(pos);
                result.rc_kmer_values.push_back(kmer);
            }
        }
    }

    // 6. Resolve per-strand thresholds
    const bool use_coverscore = (config.stage1.stage1_score_type == 1);

    // Default thresholds (non-fractional mode)
    result.resolved_threshold_fwd = config.stage1.min_stage1_score;
    result.resolved_threshold_rc = config.stage1.min_stage1_score;

    if (config.min_stage1_score_frac > 0) {
        // Fractional threshold resolution

        // Count Nqkmer per strand: count distinct positions (handles degenerate expansion)
        auto count_nqkmer = [&](const std::vector<std::pair<uint32_t, KmerInt>>& kmers) -> uint32_t {
            std::unordered_set<uint32_t> positions;
            for (const auto& [pos, kmer] : kmers) positions.insert(pos);
            return static_cast<uint32_t>(positions.size());
        };

        uint32_t Nqkmer_fwd = count_nqkmer(fwd_kmers);
        uint32_t Nqkmer_rc = count_nqkmer(rc_kmers);

        // Count Nhighfreq per strand: count position only if ALL expanded
        // k-mers at that position are high-freq (adjacency guaranteed by extract_kmers)
        auto count_nhighfreq = [&](const std::vector<std::pair<uint32_t, KmerInt>>& kmers) -> uint32_t {
            uint32_t count = 0;
            size_t i = 0;
            while (i < kmers.size()) {
                uint32_t cur_pos = kmers[i].first;
                bool all_highfreq = true;
                while (i < kmers.size() && kmers[i].first == cur_pos) {
                    if (highfreq_set.count(static_cast<uint32_t>(kmers[i].second)) == 0)
                        all_highfreq = false;
                    i++;
                }
                if (all_highfreq) count++;
            }
            return count;
        };

        uint32_t Nhighfreq_fwd = count_nhighfreq(fwd_kmers);
        uint32_t Nhighfreq_rc = count_nhighfreq(rc_kmers);

        // Compute threshold per strand
        auto resolve_threshold = [&](uint32_t Nqkmer, uint32_t Nhighfreq,
                                     const char* strand_name) -> uint32_t {
            int32_t threshold = static_cast<int32_t>(
                std::ceil(static_cast<double>(Nqkmer) * config.min_stage1_score_frac))
                - static_cast<int32_t>(Nhighfreq);

            if (threshold <= 0) {
                std::fprintf(stderr,
                    "Warning: fractional min_stage1_score threshold <= 0 "
                    "(strand=%s, Nqkmer=%u, Nhighfreq=%u, P=%.4f)\n",
                    strand_name, Nqkmer, Nhighfreq,
                    config.min_stage1_score_frac);
                return 0; // signals: skip this strand
            }
            return static_cast<uint32_t>(threshold);
        };

        result.resolved_threshold_fwd = resolve_threshold(Nqkmer_fwd, Nhighfreq_fwd, "fwd");
        result.resolved_threshold_rc = resolve_threshold(Nqkmer_rc, Nhighfreq_rc, "rc");
    }

    // 7. Resolve effective_min_score per strand
    if (config.stage2.min_score > 0) {
        // Explicit min_score: use directly for both strands
        result.effective_min_score_fwd = config.stage2.min_score;
        result.effective_min_score_rc = config.stage2.min_score;
    } else {
        // Adaptive: min_score=0 -> use resolved Stage 1 threshold
        if (config.min_stage1_score_frac > 0) {
            // Fractional mode: per-query adaptive threshold
            result.effective_min_score_fwd = result.resolved_threshold_fwd;
            result.effective_min_score_rc = result.resolved_threshold_rc;
        } else {
            // Absolute mode: use the configured min_stage1_score
            result.effective_min_score_fwd = config.stage1.min_stage1_score;
            result.effective_min_score_rc = config.stage1.min_stage1_score;
        }
    }

    return result;
}

// Explicit template instantiations
template QueryKmerData<uint16_t> preprocess_query<uint16_t>(
    const std::string&, int,
    const std::vector<const KixReader*>&,
    const KhxReader*,
    const SearchConfig&,
    uint8_t,
    const std::vector<uint32_t>&,
    const std::vector<uint16_t>&);
template QueryKmerData<uint32_t> preprocess_query<uint32_t>(
    const std::string&, int,
    const std::vector<const KixReader*>&,
    const KhxReader*,
    const SearchConfig&,
    uint8_t,
    const std::vector<uint32_t>&,
    const std::vector<uint32_t>&);

} // namespace ikafssn
