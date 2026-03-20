#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "core/types.hpp"
#include "search/stage1_filter.hpp"
#include "search/stage2_chaining.hpp"
#include "search/query_preprocessor.hpp"

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
    uint8_t  mode = 2;          // 1 = stage1 only, 2 = stage1+stage2, 3 = stage1+stage2+stage3
    uint8_t  sort_score = 2;    // 1 = stage1 score, 2 = chainscore, 3 = alnscore
    int8_t   strand = 2;       // 1 = plus only, -1 = minus only, 2 = both
    uint8_t  accept_qdegen = 1; // 0 = skip queries with degenerate bases, 1 = accept (expand)
    double min_stage1_score_frac = 0; // 0 = disabled, 0 < P < 1 = fractional mode
    uint16_t max_degen_expand = 16;  // max degenerate expansion per k-mer (0/1: disable)
    uint8_t  t = 0;  // template length (0 = contiguous)
};

struct SearchResult {
    std::string query_id;
    std::vector<ChainResult> hits;
};

// Search a single volume using pre-processed query k-mer data.
// High-freq k-mers have already been removed and thresholds resolved globally.
// buf: optional thread-local Stage1Buffer to avoid per-call allocation.
template <typename KmerInt>
SearchResult search_volume(
    const std::string& query_id,
    const QueryKmerData<KmerInt>& qdata,
    int k,
    const KixReader& kix,
    const KpxReader& kpx,
    const KsxReader& ksx,
    const OidFilter& filter,
    const SearchConfig& config,
    Stage1Buffer* buf = nullptr);

extern template SearchResult search_volume<uint16_t>(
    const std::string&, const QueryKmerData<uint16_t>&, int,
    const KixReader&, const KpxReader&, const KsxReader&,
    const OidFilter&, const SearchConfig&, Stage1Buffer*);
extern template SearchResult search_volume<uint32_t>(
    const std::string&, const QueryKmerData<uint32_t>&, int,
    const KixReader&, const KpxReader&, const KsxReader&,
    const OidFilter&, const SearchConfig&, Stage1Buffer*);

// Search a single volume using merged coding+optimal indexes ("both" mode).
// Two separate QueryKmerData are provided: one for coding, one for optimal.
// Stage 1 scores are summed across both templates; Stage 2 hits are merged.
template <typename KmerInt>
SearchResult search_volume_both(
    const std::string& query_id,
    const QueryKmerData<KmerInt>& qdata_cod,
    const QueryKmerData<KmerInt>& qdata_opt,
    int k,
    const KixReader& kix_cod, const KpxReader& kpx_cod,
    const KixReader& kix_opt, const KpxReader& kpx_opt,
    const KsxReader& ksx,
    const OidFilter& filter,
    const SearchConfig& config,
    Stage1Buffer* buf_cod = nullptr,
    Stage1Buffer* buf_opt = nullptr);

extern template SearchResult search_volume_both<uint16_t>(
    const std::string&,
    const QueryKmerData<uint16_t>&, const QueryKmerData<uint16_t>&, int,
    const KixReader&, const KpxReader&,
    const KixReader&, const KpxReader&,
    const KsxReader&, const OidFilter&, const SearchConfig&,
    Stage1Buffer*, Stage1Buffer*);
extern template SearchResult search_volume_both<uint32_t>(
    const std::string&,
    const QueryKmerData<uint32_t>&, const QueryKmerData<uint32_t>&, int,
    const KixReader&, const KpxReader&,
    const KixReader&, const KpxReader&,
    const KsxReader&, const OidFilter&, const SearchConfig&,
    Stage1Buffer*, Stage1Buffer*);

} // namespace ikafssn
