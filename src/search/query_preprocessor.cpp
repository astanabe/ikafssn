#include "search/query_preprocessor.hpp"
#include "search/volume_searcher.hpp"
#include "index/khx_reader.hpp"
#include "core/kmer_encoding.hpp"
#include "core/kmer_revcomp_simd.hpp"
#include "core/spaced_seed.hpp"
#include "core/config.hpp"

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
        max_expansion);
    return kmers;
}

template <typename KmerInt>
QueryKmerData<KmerInt> preprocess_query(
    const std::string& query_seq, int k,
    const KhxReader* khx,
    const SearchConfig& config,
    uint8_t t,
    const std::vector<uint32_t>& masks) {

    QueryKmerData<KmerInt> result;

    // 1. Extract forward k-mers
    std::vector<std::pair<uint32_t, KmerInt>> fwd_kmers;
    if (t > 0 && !masks.empty()) {
        fwd_kmers = extract_kmers_spaced<KmerInt>(query_seq, k, masks,
                        static_cast<int>(t), &result.has_multi_degen,
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
                          static_cast<int>(t), &result.has_multi_degen,
                          static_cast<int>(config.max_degen_expand));
        int span = static_cast<int>(t);
        rc_kmers.reserve(rc_raw.size());
        for (const auto& [pos, kmer] : rc_raw) {
            // Remap position: rc position p corresponds to fwd position (len - p - span)
            uint32_t fwd_pos = static_cast<uint32_t>(query_seq.size()) - pos - static_cast<uint32_t>(span);
            rc_kmers.emplace_back(fwd_pos, kmer);
        }
    } else {
        // Contiguous: SIMD batch revcomp via SoA staging buffers.
        // The destination is a vector<pair> (AoS), but extracting kmer values
        // into a contiguous buffer lets the SIMD kernel run on aligned dense
        // data; the AoS reassembly is a cheap two-vector zip.
        const std::size_t nfwd = fwd_kmers.size();
        std::vector<KmerInt> tmp_in(nfwd), tmp_out(nfwd);
        for (std::size_t i = 0; i < nfwd; i++) tmp_in[i] = fwd_kmers[i].second;
        kmer_revcomp_batch<KmerInt>(tmp_in.data(), tmp_out.data(), nfwd, k);
        rc_kmers.resize(nfwd);
        for (std::size_t i = 0; i < nfwd; i++) {
            rc_kmers[i] = {fwd_kmers[i].first, tmp_out[i]};
        }
    }

    // 3. SoA conversion (no search-time high-freq filtering; build-time
    //    exclusion via .khx is the only frequency gate now).
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
    to_soa(fwd_kmers, result.fwd_positions, result.fwd_kmer_values);
    to_soa(rc_kmers, result.rc_positions, result.rc_kmer_values);

    // 4. Resolve per-strand thresholds
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
        // k-mers at that position are excluded by .khx (build-time exclusion).
        // Adjacency guaranteed by extract_kmers.
        auto count_nhighfreq = [&](const std::vector<std::pair<uint32_t, KmerInt>>& kmers) -> uint32_t {
            if (khx == nullptr) return 0;
            uint32_t count = 0;
            size_t i = 0;
            while (i < kmers.size()) {
                uint32_t cur_pos = kmers[i].first;
                bool all_excluded = true;
                while (i < kmers.size() && kmers[i].first == cur_pos) {
                    if (!khx->is_excluded(static_cast<uint32_t>(kmers[i].second)))
                        all_excluded = false;
                    i++;
                }
                if (all_excluded) count++;
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

    // 5. Resolve effective_min_score per strand
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
    const KhxReader*,
    const SearchConfig&,
    uint8_t,
    const std::vector<uint32_t>&);
template QueryKmerData<uint32_t> preprocess_query<uint32_t>(
    const std::string&, int,
    const KhxReader*,
    const SearchConfig&,
    uint8_t,
    const std::vector<uint32_t>&);

} // namespace ikafssn
