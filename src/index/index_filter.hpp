#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace ikafssn {

class Logger;

// Cross-volume frequency filtering for -max_freq_build.
//
// Reads .kix.tmp files from all volumes, sums per-kmer counts across volumes,
// determines which k-mers exceed freq_threshold, then:
//   - Writes filtered .kix/.kpx from each .kix.tmp/.kpx.tmp (excluded k-mers removed)
//   - Renames .ksx.tmp -> .ksx
//   - Generates shared .khx at khx_path
//   - Removes .tmp files on success
//
// vol_prefixes: per-volume output prefixes (e.g. "/out/nt.00.11mer", "/out/nt.01.11mer")
// khx_path: shared .khx file path (e.g. "/out/nt.11mer.khx")
// k: k-mer length
// freq_threshold: exclude k-mers with global count > freq_threshold
//
// Returns true on success.
bool filter_volumes_cross_volume(
    const std::vector<std::string>& vol_prefixes,
    const std::string& khx_path,
    int k,
    uint64_t freq_threshold,
    int filter_threads,
    const Logger& logger);

// Phase 7d: structural validation for a (.kix, .kpx) volume pair.
//
// Walks every k-mer's posting list and verifies that the byte length
// recorded in the .kpx EF dictionary matches the bytes actually
// consumed by the kind map + partition groups + short1 FOR stream +
// short2 occ_count[] + short2 FOR stream.  Catches silent kind-map
// or FOR-stream corruption that the v9 dedup'd headers can no longer
// detect by redundancy.  Also returns the recomputed total position
// count via `total_position_count_out` so the caller can populate the
// .kpx header for filtered builds (the legacy v8 recompute path read
// the wrong field).
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
