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
#include "search/stage3_alignment.hpp"
#include "protocol/messages.hpp"

#include <tbb/task_arena.h>

namespace ikafssn {

class Server;  // forward declaration
struct DatabaseEntry;  // forward declaration

// Pre-opened volume data (shared read-only across threads)
struct ServerVolumeData {
    KixReader kix;
    KpxReader kpx;
    KsxReader ksx;
    uint16_t volume_index;
    uint64_t total_bases = 0;
};

// A group of volumes for a specific k-mer size
struct KmerGroup {
    int k;
    uint8_t kmer_type; // 0 = uint16_t, 1 = uint32_t
    std::vector<ServerVolumeData> volumes;
    KhxReader khx;  // shared .khx for this k-mer size
};

// Process a search request using loaded index data from a specific database.
// Acquires per-sequence permits via server semaphore; rejected queries
// are returned in resp.rejected_query_ids for client retry.
SearchResponse process_search_request(
    const SearchRequest& req,
    const DatabaseEntry& db,
    Server& server,
    tbb::task_arena& arena);

} // namespace ikafssn
