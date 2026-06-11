#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace ikafssn {

class Logger;

// Write a .khx file from a pre-computed exclusion bitset.
// excluded[i] == true means k-mer i is excluded.
// Returns true on success.
bool write_khx_bitset(const std::string& path, int k,
                      const std::vector<bool>& excluded,
                      const Logger& logger);

} // namespace ikafssn
