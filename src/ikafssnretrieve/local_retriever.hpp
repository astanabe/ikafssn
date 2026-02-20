#pragma once

#include <cstdint>
#include <ostream>
#include <string>
#include <vector>

#include "io/result_writer.hpp"

namespace ikafssn {

class BlastDbReader;

struct RetrieveOptions {
    uint32_t context = 0;  // bases to add before/after match region
};

// Retrieve matched subsequences from a local BLAST DB.
// Writes FASTA records to the output stream.
// Returns the number of successfully retrieved sequences.
uint32_t retrieve_local(const std::vector<OutputHit>& hits,
                        const std::string& db_path,
                        const RetrieveOptions& opts,
                        std::ostream& out);

} // namespace ikafssn
