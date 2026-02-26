#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "io/result_writer.hpp"
#include "io/fasta_reader.hpp"
#include "util/logger.hpp"

namespace ikafssn {

struct Stage3Config {
    int gapopen = 10;
    int gapext = 1;
    bool traceback = false;
    double min_pident = 0.0;
    uint32_t min_nident = 0;
    int fetch_threads = 8;   // threads for BLAST DB sequence fetch
};

// Run Stage 3 alignment on merged OutputHits.
// - hits: Stage 2 results (modified in-place with alignment data)
// - queries: original FASTA query sequences
// - db_path: BLAST DB path for subject sequence retrieval
// - context_is_ratio/context_ratio/context_abs: -context option values
// - Fetch thread count is controlled by config.fetch_threads
// Returns filtered hits (min_pident/min_nident applied).
std::vector<OutputHit> run_stage3(
    std::vector<OutputHit>& hits,
    const std::vector<FastaRecord>& queries,
    const std::string& db_path,
    const Stage3Config& config,
    bool context_is_ratio,
    double context_ratio,
    uint32_t context_abs,
    const Logger& logger);

} // namespace ikafssn
