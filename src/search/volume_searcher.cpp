#include "search/volume_searcher.hpp"
#include "search/oid_filter.hpp"
#include "search/seq_id_decoder.hpp"
#include "search/posting_decoder.hpp"
#include "search/stage1_filter.hpp"
#include "search/stage2_chaining.hpp"
#include "index/kix_reader.hpp"
#include "index/kpx_reader.hpp"
#include "index/ksx_reader.hpp"
#include "core/kmer_encoding.hpp"

#include <algorithm>
#include <unordered_set>

namespace ikafssn {

template <typename KmerInt>
static std::vector<std::pair<uint32_t, KmerInt>>
extract_kmers(const std::string& seq, int k) {
    std::vector<std::pair<uint32_t, KmerInt>> kmers;
    KmerScanner<KmerInt> scanner(k);
    scanner.scan(seq.data(), seq.size(), [&](uint32_t pos, KmerInt kmer) {
        kmers.emplace_back(pos, kmer);
    });
    return kmers;
}

template <typename KmerInt>
static std::vector<ChainResult>
search_one_strand(const std::vector<std::pair<uint32_t, KmerInt>>& query_kmers,
                  int k,
                  bool is_reverse,
                  const KixReader& kix,
                  const KpxReader& kpx,
                  const OidFilter& filter,
                  const SearchConfig& config) {

    // Stage 1: candidate selection
    std::vector<SeqId> candidates = stage1_filter(query_kmers, kix, filter, config.stage1);
    if (candidates.empty()) return {};

    // Build candidate set for fast lookup
    std::unordered_set<SeqId> candidate_set(candidates.begin(), candidates.end());

    // Compute effective max_freq (same logic as stage1)
    uint32_t max_freq = config.stage1.max_freq;
    if (max_freq == 0) {
        double mean = static_cast<double>(kix.total_postings()) /
                      static_cast<double>(kix.table_size());
        max_freq = static_cast<uint32_t>(mean * 10.0);
        if (max_freq < 1000) max_freq = 1000;
        if (max_freq > 100000) max_freq = 100000;
    }

    // Stage 2: collect hits for candidates
    std::unordered_map<SeqId, std::vector<Hit>> hits_per_seq;

    const uint64_t* offsets = kix.offsets();
    const uint32_t* counts = kix.counts();
    const uint8_t* id_data = kix.posting_data();
    const uint64_t* pos_offsets = kpx.pos_offsets();
    const uint8_t* pos_data = kpx.posting_data();

    for (const auto& [q_pos, kmer] : query_kmers) {
        uint64_t kmer_idx = static_cast<uint64_t>(kmer);
        uint32_t cnt = counts[kmer_idx];
        if (cnt == 0 || cnt > max_freq) continue;

        SeqIdDecoder id_decoder(id_data + offsets[kmer_idx]);
        PosDecoder pos_decoder(pos_data + pos_offsets[kmer_idx]);

        for (uint32_t i = 0; i < cnt; i++) {
            SeqId sid = id_decoder.next();
            uint32_t s_pos = pos_decoder.next(id_decoder.was_new_seq());

            if (candidate_set.count(sid)) {
                hits_per_seq[sid].push_back({q_pos, s_pos});
            }
        }
    }

    // Chain hits for each candidate
    std::vector<ChainResult> results;
    for (SeqId sid : candidates) {
        auto it = hits_per_seq.find(sid);
        if (it == hits_per_seq.end()) continue;

        ChainResult cr = chain_hits(it->second, sid, k, is_reverse, config.stage2);
        if (cr.score >= config.stage2.min_score) {
            results.push_back(cr);
        }
    }

    return results;
}

template <typename KmerInt>
SearchResult search_volume(
    const std::string& query_id,
    const std::string& query_seq,
    int k,
    const KixReader& kix,
    const KpxReader& kpx,
    const KsxReader& ksx,
    const OidFilter& filter,
    const SearchConfig& config) {

    SearchResult result;
    result.query_id = query_id;

    // Extract k-mers from forward strand
    auto fwd_kmers = extract_kmers<KmerInt>(query_seq, k);

    // Search forward strand
    auto fwd_results = search_one_strand(fwd_kmers, k, false, kix, kpx, filter, config);
    result.hits.insert(result.hits.end(), fwd_results.begin(), fwd_results.end());

    // Generate reverse complement k-mers
    std::vector<std::pair<uint32_t, KmerInt>> rc_kmers;
    rc_kmers.reserve(fwd_kmers.size());
    for (const auto& [pos, kmer] : fwd_kmers) {
        KmerInt rc = kmer_revcomp(kmer, k);
        // For reverse complement search, the query position convention:
        // q_pos in the reverse strand maps to the same position in forward coord
        rc_kmers.emplace_back(pos, rc);
    }

    // Search reverse complement
    auto rc_results = search_one_strand(rc_kmers, k, true, kix, kpx, filter, config);
    result.hits.insert(result.hits.end(), rc_results.begin(), rc_results.end());

    // Sort all hits by score descending
    std::sort(result.hits.begin(), result.hits.end(),
              [](const ChainResult& a, const ChainResult& b) {
                  return a.score > b.score;
              });

    // Truncate to num_results
    if (result.hits.size() > config.num_results) {
        result.hits.resize(config.num_results);
    }

    return result;
}

template SearchResult search_volume<uint16_t>(
    const std::string&, const std::string&, int,
    const KixReader&, const KpxReader&, const KsxReader&,
    const OidFilter&, const SearchConfig&);
template SearchResult search_volume<uint32_t>(
    const std::string&, const std::string&, int,
    const KixReader&, const KpxReader&, const KsxReader&,
    const OidFilter&, const SearchConfig&);

} // namespace ikafssn
