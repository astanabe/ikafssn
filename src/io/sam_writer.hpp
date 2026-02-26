#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "io/result_writer.hpp"

namespace ikafssn {

// Write results in SAM format.
// output_path: file path, or "-" for stdout.
void write_results_sam(const std::string& output_path,
                       const std::vector<OutputHit>& hits,
                       uint8_t stage1_score_type);

// Write results in BAM format.
// output_path: file path (must not be empty or "-").
void write_results_bam(const std::string& output_path,
                       const std::vector<OutputHit>& hits,
                       uint8_t stage1_score_type);

} // namespace ikafssn
