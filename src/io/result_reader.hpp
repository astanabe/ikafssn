#pragma once

#include <string>
#include <vector>

#include "io/result_writer.hpp"

namespace ikafssn {

// Parse search results from TSV (tab-delimited) format.
// Handles both file input and stdin ("-").
// Skips header lines (starting with '#') and blank lines.
// Returns empty vector on error.
std::vector<OutputHit> read_results_tsv(const std::string& path);

// Parse search results from an input stream.
std::vector<OutputHit> read_results_tsv(std::istream& in);

} // namespace ikafssn
