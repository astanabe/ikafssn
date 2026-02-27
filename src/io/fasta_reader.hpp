#pragma once

#include <cstdint>
#include <istream>
#include <string>
#include <vector>

namespace ikafssn {

struct FastaRecord {
    std::string id;       // sequence ID (first word after '>')
    std::string sequence; // concatenated sequence lines (uppercase)
};

// Read all records from an input stream.
std::vector<FastaRecord> read_fasta_stream(std::istream& in);

// Read all records from a FASTA file.
// path can be "-" for stdin.
// Returns empty vector on error.
std::vector<FastaRecord> read_fasta(const std::string& path);

} // namespace ikafssn
