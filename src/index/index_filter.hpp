#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace ikafssn {

class Logger;

// Cross-volume frequency filtering for -max_freq_build.
//
// Reads .kix.tmp files from all volumes, sums per-kmer counts across
// volumes, identifies k-mers above freq_threshold, then writes the
// filtered .kix / .kpx and .ksx at vol_prefixes_final[vi] from the
// .kix.tmp / .kpx.tmp / .ksx.tmp at vol_prefixes_tmp[vi], generates the
// shared .khx at khx_path, and removes the .tmp files on success.  The
// two prefix vectors differ when the final file name encodes parameters
// not yet known when the .tmp files were written.
bool filter_volumes_cross_volume(
    const std::vector<std::string>& vol_prefixes_tmp,
    const std::vector<std::string>& vol_prefixes_final,
    const std::string& khx_path,
    int k,
    uint64_t freq_threshold,
    int filter_threads,
    const Logger& logger);

// Structural validation for a (.kix, .kpx) volume pair.
//
// Walks every k-mer's posting list and verifies that the byte length
// recorded in the .kpx EF dictionary matches the bytes actually
// consumed by the kind map + partition groups + short1 FOR stream +
// short2 occ_count[] + short2 FOR stream.  Catches silent kind-map
// or FOR-stream corruption.  Also returns the recomputed total
// position count via `total_position_count_out` so the caller can
// populate the .kpx header for filtered builds.
//
//   kix_path / kpx_path   absolute paths to the volume's index files
//   total_position_count_out  optional; receives sum of partition
//                             occ_counts + short1_count + sum(short2
//                             occ_count[]) across the whole volume
//   logger                error / info sink
//
// Returns true if every k-mer's posting list is internally consistent.
// On failure, logs the offending k-mer index and the mismatch reason
// before returning false.  When kpx_path is empty, only the .kix side
// is validated (mode-1 indexes).
bool validate_volume(const std::string& kix_path,
                     const std::string& kpx_path,
                     uint64_t* total_position_count_out,
                     const Logger& logger);

} // namespace ikafssn
