#pragma once

#include <cstdint>
#include <ostream>
#include <string>
#include <vector>

#include "io/result_writer.hpp"

namespace ikafssn {

class BlastDbReader;

struct RetrieveOptions {
    // Context bases added before/after the match region.  In ratio mode the
    // per-hit context is derived from each hit's query length; otherwise the
    // fixed `context` value is used for every hit.
    uint32_t context = 0;        // absolute context (ratio == false)
    bool     is_ratio = false;
    double   ratio = 0.0;        // per-hit ctx = round(hit.qlen * ratio)
};

// Retrieve matched subsequences from a local BLAST DB.
// Writes FASTA records to the output stream.
// Returns the number of successfully retrieved sequences.
uint32_t retrieve_local(const std::vector<OutputHit>& hits,
                        const std::string& db_path,
                        const RetrieveOptions& opts,
                        std::ostream& out);

} // namespace ikafssn
