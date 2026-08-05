#pragma once

#include <string>
#include <vector>

#include "io/result_writer.hpp"

namespace ikafssn {

// Parse search results from TSV (tab-delimited) format.
// Handles both file input and stdin ("-").
// Skips header lines (starting with '#') and blank lines.
// Returns empty vector on error.
//
// A missing column reads back as 0, so `has_volume_column` reports whether
// the input carried one — the only way to tell it from a genuine volume 0.
// It stays false when nothing could be parsed.
std::vector<OutputHit> read_results_tsv(const std::string& path,
                                        bool* has_volume_column = nullptr);

// Parse search results from an input stream.
std::vector<OutputHit> read_results_tsv(std::istream& in,
                                        bool* has_volume_column = nullptr);

} // namespace ikafssn
