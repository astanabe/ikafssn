#pragma once

#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include "core/types.hpp"
#include "index/ksx_reader.hpp"
#include "index/khx_reader.hpp"
#include "io/volume_discovery.hpp"
#include "search/oid_filter.hpp"
#include "search/search_config.hpp"
#include "search/stage3_alignment.hpp"
#include "protocol/messages.hpp"

#include <tbb/task_arena.h>

namespace ikafssn {

class Server;  // forward declaration
struct DatabaseEntry;  // forward declaration

// Per-volume static metadata for one (k, t, template_type) group.  The
// .kix / .kpx readers are opened/closed per-request inside the search
// orchestrator (avoids concurrent-request madvise contention on shared
// readers); only .ksx is held open here for accession / seq_length lookups
// and OidFilter construction.
struct ServerVolumeData {
    DiscoveredVolume files;
    KsxReader ksx;
    uint32_t num_sequences = 0;
    uint64_t total_distinct_postings = 0;
    std::string db_name;          // KixHeader::db captured at validation
    uint16_t volume_index = 0;
    uint64_t total_bases = 0;
};

// A group of volumes for one full index variant.  Every volume in the group
// shares all 8 identifying parameters (k / t / template_type and the five
// indexing parameters below); different variants live in different groups.
struct KmerGroup {
    int k;
    uint8_t kmer_type; // 0 = uint16_t, 1 = uint32_t
    uint8_t t = 0;              // template length
    uint8_t template_type = 0;  // TemplateType enum value
    uint32_t min_seq_length = 0;
    uint32_t min_length_split = 0;
    uint32_t overlap_length = 0;
    uint64_t max_freq_build = 1;
    uint32_t max_degen_expand = 0;
    std::vector<ServerVolumeData> volumes;
    KhxReader khx;  // shared .khx for this k-mer size
};

// Process a search request using loaded index data from a specific database.
// Acquires per-sequence permits via server semaphore; rejected queries
// are returned in resp.rejected_qseqids for client retry.
SearchResponse process_search_request(
    const SearchRequest& req,
    const DatabaseEntry& db,
    Server& server,
    tbb::task_arena& arena);

} // namespace ikafssn
