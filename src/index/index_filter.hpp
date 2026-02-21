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
    const Logger& logger);

} // namespace ikafssn
