#pragma once

#include <string>
#include <vector>

namespace ikafssn {

// Read a seqidlist file (accession list for filtering).
// Supports:
//   - Text format: one accession per line, '>' prefix trimmed, blank lines skipped
//   - Binary format: blastdb_aliastool -seqid_file_in generated files
// Format is auto-detected via magic bytes.
// Returns list of accession strings.
std::vector<std::string> read_seqidlist(const std::string& path);

} // namespace ikafssn
