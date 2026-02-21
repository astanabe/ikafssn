#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace ikafssn {

class Logger;

// Write a .khx file recording which k-mers were excluded during index build.
// counts: per-kmer count array (before zeroing excluded entries).
// freq_threshold: counts[i] > freq_threshold means k-mer i is excluded.
// Returns true on success.
bool write_khx(const std::string& path, int k,
               const std::vector<uint32_t>& counts,
               uint64_t freq_threshold,
               const Logger& logger);

// Write a .khx file from a pre-computed exclusion bitset.
// excluded[i] == true means k-mer i is excluded.
// Returns true on success.
bool write_khx_bitset(const std::string& path, int k,
                      const std::vector<bool>& excluded,
                      const Logger& logger);

} // namespace ikafssn
