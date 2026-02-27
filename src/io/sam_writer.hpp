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

// Write results in the appropriate format (tab/json/sam/bam), dispatching
// to the correct writer. Returns true on success.
bool write_all_results(const std::string& output_path,
                       const std::vector<OutputHit>& hits,
                       OutputFormat fmt,
                       uint8_t mode,
                       uint8_t stage1_score_type,
                       bool stage3_traceback);

// Merge multiple SAM batch files into a single output.
// @SQ headers are unioned across all files; @HD and @PG from the first file.
// Records are remapped to the merged header's tid space.
// as_bam=true writes BAM ("wb"), false writes SAM ("w").
// Returns true on success.
bool merge_sam_files(const std::vector<std::string>& batch_paths,
                     const std::string& output_path, bool as_bam);

} // namespace ikafssn
