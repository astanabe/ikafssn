#pragma once

#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include "core/types.hpp"
#include "index/kix_reader.hpp"
#include "index/kpx_reader.hpp"
#include "index/ksx_reader.hpp"
#include "index/khx_reader.hpp"
#include "search/oid_filter.hpp"
#include "search/volume_searcher.hpp"
#include "protocol/messages.hpp"

namespace ikafssn {

// Pre-opened volume data (shared read-only across threads)
struct ServerVolumeData {
    KixReader kix;
    KpxReader kpx;
    KsxReader ksx;
    KhxReader khx;
    uint16_t volume_index;
};

// A group of volumes for a specific k-mer size
struct KmerGroup {
    int k;
    uint8_t kmer_type; // 0 = uint16_t, 1 = uint32_t
    std::vector<ServerVolumeData> volumes;
};

// Process a search request using loaded index data.
// Returns a SearchResponse.
SearchResponse process_search_request(
    const SearchRequest& req,
    const std::map<int, KmerGroup>& kmer_groups,
    int default_k,
    const SearchConfig& default_config);

} // namespace ikafssn
