#pragma once

// Single-volume search convenience wrappers used only by the test suite.
// Production search runs through search_orchestrator::run_search(); these
// helpers expose the per-volume Stage 1 / Stage 2 building blocks as one
// SearchResult so individual stages can be exercised in isolation.

#include <cstdint>
#include <string>
#include <vector>

#include "core/types.hpp"
#include "search/search_config.hpp"
#include "search/stage1_filter.hpp"
#include "search/stage2_chaining.hpp"
#include "search/query_preprocessor.hpp"

namespace ikafssn {

class KixReader;
class KpxReader;
class KsxReader;
class KhxReader;
class OidFilter;

struct SearchResult {
    std::string query_id;
    std::vector<ChainResult> hits;
};

// Search a single volume using pre-processed query k-mer data.
// High-freq k-mers have already been removed and thresholds resolved globally.
// buf: thread-local Stage1Buffer (must have buf.width set) to avoid per-call
// allocation. Required.
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
    Stage1Buffer& buf);

extern template SearchResult search_volume<uint16_t>(
    const std::string&, const QueryKmerData<uint16_t>&, int,
    const KixReader&, const KpxReader&, const KsxReader&,
    const OidFilter&, const SearchConfig&, Stage1Buffer&);
extern template SearchResult search_volume<uint32_t>(
    const std::string&, const QueryKmerData<uint32_t>&, int,
    const KixReader&, const KpxReader&, const KsxReader&,
    const OidFilter&, const SearchConfig&, Stage1Buffer&);

// Search a single volume using merged coding+optimal indexes ("both" mode).
// Two separate QueryKmerData are provided: one for coding, one for optimal.
// Stage 1 runs both templates against the *same* Stage1Buffer so the
// per-(sid, q_pos) dedup carries across templates; Stage 2 hits are merged
// across both indexes before chaining.
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
    Stage1Buffer& buf);

extern template SearchResult search_volume_both<uint16_t>(
    const std::string&,
    const QueryKmerData<uint16_t>&, const QueryKmerData<uint16_t>&, int,
    const KixReader&, const KpxReader&,
    const KixReader&, const KpxReader&,
    const KsxReader&, const OidFilter&, const SearchConfig&,
    Stage1Buffer&);
extern template SearchResult search_volume_both<uint32_t>(
    const std::string&,
    const QueryKmerData<uint32_t>&, const QueryKmerData<uint32_t>&, int,
    const KixReader&, const KpxReader&,
    const KixReader&, const KpxReader&,
    const KsxReader&, const OidFilter&, const SearchConfig&,
    Stage1Buffer&);

} // namespace ikafssn
