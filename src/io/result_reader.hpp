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
// A missing column reads back as 0, so callers that need to tell "volume 0"
// from "no volume column at all" pass `has_volume_column`; it is set to
// whether the input actually carried the column.  It stays false when
// nothing could be parsed.
std::vector<OutputHit> read_results_tsv(const std::string& path,
                                        bool* has_volume_column = nullptr);

// Parse search results from an input stream.
std::vector<OutputHit> read_results_tsv(std::istream& in,
                                        bool* has_volume_column = nullptr);

} // namespace ikafssn
