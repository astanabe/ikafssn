#include "search/query_preprocessor.hpp"
#include "search/volume_searcher.hpp"
#include "index/khx_reader.hpp"
#include "core/kmer_encoding.hpp"
#include "core/kmer_revcomp_simd.hpp"
#include "core/spaced_seed.hpp"
#include "core/config.hpp"
#include "protocol/messages.hpp"

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

// Detect first non-IUPAC, non-ATGC character in the sequence.
// Returns string::npos if none, otherwise the first position.
static size_t find_invalid_char(const std::string& seq) {
    for (size_t i = 0; i < seq.size(); i++) {
        unsigned char c = static_cast<unsigned char>(seq[i]);
        uint8_t enc = encode_base(static_cast<char>(c));
        uint8_t ncbi4na = degenerate_ncbi4na_table()[c];
        if (enc == BASE_ENCODE_INVALID && ncbi4na == 0) return i;
    }
    return std::string::npos;
}

template <typename KmerInt>
QueryKmerData<KmerInt> preprocess_query(
    const std::string& query_seq, int k,
    const KhxReader* khx,
    const SearchConfig& config,
    uint8_t t,
    const std::vector<uint32_t>& masks) {

    QueryKmerData<KmerInt> result;
    result.qlen = static_cast<uint32_t>(query_seq.size());

    // 0. Length check.
    // First, the explicit -min_query_length floor (v10).  Queries shorter
    // than this are skipped before the span-based check below.
    if (config.min_query_length > 0 &&
        query_seq.size() < config.min_query_length) {
        result.skip_reason = kSkipQueryTooShort;
        char buf[128];
        std::snprintf(buf, sizeof(buf), "length=%zu < min_query_length=%u",
                      query_seq.size(),
                      static_cast<unsigned>(config.min_query_length));
        result.skip_detail = buf;
        return result;
    }
    // Then the v10 (Phase 3) overlap_length ceiling.  The index reports its
    // overlap_length via the .ksx header; the orchestrator copies it into
    // SearchConfig::max_query_length.  Queries longer than the overlap are
    // skipped because the Stage 2/3 dedup keys (parent-relative coordinates)
    // assume every chain hit fits inside at most two adjacent fragments.
    // max_query_length == 0 disables the check (degenerate fragment-table
    // case: min_length_split == 0, no splitting).
    if (config.max_query_length > 0 &&
        query_seq.size() > config.max_query_length) {
        result.skip_reason = kSkipQueryTooLong;
        char buf[160];
        std::snprintf(buf, sizeof(buf),
                      "length=%zu > overlap_length=%u; queries longer than the "
                      "index's overlap_length are not supported",
                      query_seq.size(),
                      static_cast<unsigned>(config.max_query_length));
        result.skip_detail = buf;
        return result;
    }
    // The window-count Nqkmer also requires seq_len >= span.
    const int span = (t > 0) ? static_cast<int>(t) : k;
    if (static_cast<int>(query_seq.size()) < span) {
        result.skip_reason = kSkipQueryTooShort;
        char buf[128];
        std::snprintf(buf, sizeof(buf), "length=%zu < %s=%d",
                      query_seq.size(), (t > 0) ? "t" : "k", span);
        result.skip_detail = buf;
        return result;
    }

    // 0a. accept_qdegen=0 rejection (centralized here)
    if (config.accept_qdegen == 0 && contains_degenerate_base(query_seq)) {
        result.skip_reason = kSkipDegenRejected;
        result.skip_detail = "query contains IUPAC degenerate bases";
        return result;
    }

    // 0b. Truly invalid character detection (e.g. '*', '!')
    {
        size_t bad = find_invalid_char(query_seq);
        if (bad != std::string::npos) {
            result.skip_reason = kSkipInvalidChar;
            char buf[160];
            std::snprintf(buf, sizeof(buf), "first invalid char '%c' at pos %zu",
                          query_seq[bad], bad);
            result.skip_detail = buf;
            return result;
        }
    }

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

    // 2. Build reverse complement k-mers
    std::vector<std::pair<uint32_t, KmerInt>> rc_kmers;
    if (t > 0 && !masks.empty()) {
        // Spaced seed: scan RC string with same templates, remap positions
        std::string rc_seq = reverse_complement_string(query_seq);
        auto rc_raw = extract_kmers_spaced<KmerInt>(rc_seq, k, masks,
                          static_cast<int>(t), &result.has_multi_degen,
                          static_cast<int>(config.max_degen_expand));
        rc_kmers.reserve(rc_raw.size());
        for (const auto& [pos, kmer] : rc_raw) {
            // Remap position: rc position p corresponds to fwd position (len - p - span)
            uint32_t fwd_pos = static_cast<uint32_t>(query_seq.size()) - pos - static_cast<uint32_t>(span);
            rc_kmers.emplace_back(fwd_pos, kmer);
        }
    } else if (!fwd_kmers.empty()) {
        // Contiguous: SIMD batch revcomp via SoA staging buffers.
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

    // 4. Resolve per-strand thresholds.
    //
    // Spec (Phase 1 fix):
    //   Nqkmer    = max(0, seq_len - span + 1)   (pure window count)
    //   Nhighfreq = #{p : ANY emitted k-mer at p is .khx-excluded}     (case 1)
    //             + (Nqkmer - #emitted positions)                      (case 2 + case 3)
    //   threshold = ceil(Nqkmer * P) - Nhighfreq
    //   threshold < 1 => kSkipThresholdUnreachable
    const uint32_t Nqkmer = static_cast<uint32_t>(query_seq.size() - span + 1);

    // Default thresholds (non-fractional mode)
    result.resolved_threshold_fwd = config.stage1.min_stage1_score;
    result.resolved_threshold_rc = config.stage1.min_stage1_score;

    if (config.min_stage1_score_frac > 0) {
        // Count distinct emitted positions and ANY-excluded positions for each
        // strand. The rc remap leaves the per-position grouping intact (same
        // positions, just relabeled to fwd-coords), so the same logic works.
        auto compute_nhighfreq = [&](const std::vector<std::pair<uint32_t, KmerInt>>& kmers) -> uint32_t {
            // Sort-merge by position to detect ANY-excluded per cluster.
            // Inputs from extract_kmers_* are already grouped by emission order
            // (per-position cluster), so we can scan linearly.
            std::unordered_set<uint32_t> emitted;
            uint32_t any_excluded = 0;
            size_t i = 0;
            while (i < kmers.size()) {
                uint32_t cur_pos = kmers[i].first;
                emitted.insert(cur_pos);
                bool any_ex = false;
                while (i < kmers.size() && kmers[i].first == cur_pos) {
                    if (khx != nullptr &&
                        khx->is_excluded(static_cast<uint32_t>(kmers[i].second))) {
                        any_ex = true;
                    }
                    i++;
                }
                if (any_ex) any_excluded++;
            }
            uint32_t emitted_count = static_cast<uint32_t>(emitted.size());
            uint32_t unemit = (Nqkmer >= emitted_count) ? (Nqkmer - emitted_count) : 0;
            return any_excluded + unemit;
        };

        uint32_t Nhighfreq_fwd = compute_nhighfreq(fwd_kmers);
        uint32_t Nhighfreq_rc = compute_nhighfreq(rc_kmers);

        auto resolve_threshold = [&](uint32_t Nhighfreq) -> int32_t {
            int32_t threshold = static_cast<int32_t>(
                std::ceil(static_cast<double>(Nqkmer) * config.min_stage1_score_frac))
                - static_cast<int32_t>(Nhighfreq);
            return threshold;
        };

        int32_t th_fwd = resolve_threshold(Nhighfreq_fwd);
        int32_t th_rc = resolve_threshold(Nhighfreq_rc);

        // If both strands fall below 1, the query as a whole is skipped.
        // (Single-strand searches collapse one of the two effectively.)
        bool fwd_ok = (config.strand != -1) && (th_fwd >= 1);
        bool rc_ok  = (config.strand !=  1) && (th_rc  >= 1);
        if (!fwd_ok && !rc_ok) {
            result.skip_reason = kSkipThresholdUnreachable;
            char buf[256];
            std::snprintf(buf, sizeof(buf),
                "Nqkmer=%u Nhighfreq_fwd=%u Nhighfreq_rc=%u P=%.4f "
                "threshold_fwd=%d threshold_rc=%d",
                Nqkmer, Nhighfreq_fwd, Nhighfreq_rc, config.min_stage1_score_frac,
                th_fwd, th_rc);
            result.skip_detail = buf;
            return result;
        }

        result.resolved_threshold_fwd = (th_fwd >= 1) ? static_cast<uint32_t>(th_fwd) : 0;
        result.resolved_threshold_rc  = (th_rc  >= 1) ? static_cast<uint32_t>(th_rc)  : 0;
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
