#pragma once

#include <cstdint>
#include <string>

namespace ikafssn {

class BlastDbReader;
class Logger;

// Configuration for index building.
struct IndexBuilderConfig {
    int k = 11;                         // k-mer length
    uint64_t buffer_size = uint64_t(8) << 30; // 8 GB default
    int partitions = 4;                 // number of partitions (power of 2)
    uint64_t max_freq_build = 0;        // 0 = no exclusion
    int threads = 1;                    // threads (counting + partition scan + sort)
    bool verbose = false;
};

// Build .kix, .kpx, .ksx index files for a single BLAST DB volume.
// Template parameter KmerInt selects uint16_t (k<=8) or uint32_t (k>=9).
//
// output_prefix: directory + base name, e.g. "/out/nt.00.11mer"
//   -> writes nt.00.11mer.ksx, nt.00.11mer.kix, nt.00.11mer.kpx
//
// volume_index / total_volumes: metadata stored in kix header.
//
// Returns true on success.
template <typename KmerInt>
bool build_index(BlastDbReader& db,
                 const IndexBuilderConfig& config,
                 const std::string& output_prefix,
                 uint16_t volume_index,
                 uint16_t total_volumes,
                 const std::string& db_name,
                 const Logger& logger);

} // namespace ikafssn
