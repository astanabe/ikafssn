#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "core/types.hpp"
#include "search/stage1_filter.hpp"
#include "search/stage2_chaining.hpp"

namespace ikafssn {

class KixReader;
class KpxReader;
class KsxReader;
class KhxReader;
class OidFilter;

struct SearchConfig {
    Stage1Config stage1;
    Stage2Config stage2;
    uint32_t num_results = 0;   // max results per query (0 = unlimited)
    uint8_t  mode = 2;          // 1 = stage1 only, 2 = stage1+stage2
    uint8_t  sort_score = 2;    // 1 = stage1 score, 2 = chainscore
    double min_stage1_score_frac = 0; // 0 = disabled, 0 < P < 1 = fractional mode
};

struct SearchResult {
    std::string query_id;
    std::vector<ChainResult> hits;
};

// Search a single volume for a single query (both strands).
// Template parameter KmerInt: uint16_t or uint32_t.
// When mode=1, kpx is not accessed (may be unopened).
// khx: nullable pointer to KhxReader for build-time exclusion info.
template <typename KmerInt>
SearchResult search_volume(
    const std::string& query_id,
    const std::string& query_seq,
    int k,
    const KixReader& kix,
    const KpxReader& kpx,
    const KsxReader& ksx,
    const OidFilter& filter,
    const SearchConfig& config,
    const KhxReader* khx = nullptr);

extern template SearchResult search_volume<uint16_t>(
    const std::string&, const std::string&, int,
    const KixReader&, const KpxReader&, const KsxReader&,
    const OidFilter&, const SearchConfig&, const KhxReader*);
extern template SearchResult search_volume<uint32_t>(
    const std::string&, const std::string&, int,
    const KixReader&, const KpxReader&, const KsxReader&,
    const OidFilter&, const SearchConfig&, const KhxReader*);

} // namespace ikafssn
