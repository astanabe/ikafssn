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
    double min_ppositive = 0.0;
    uint32_t min_npositive = 0;
    std::string score_matrix = "degmatch";
    // Heap budget for Stage 3 batch loop.  When 0 (default), the entire
    // hit set is processed in a single batch.  When > 0, hits are grouped
    // by (qseqid, sseqid, sstrand) and bin-packed into batches whose
    // estimated heap footprint stays under this budget.  A single
    // oversize group falls back to a solo batch.
    uint64_t posting_budget = 0;
};

// Run Stage 3 alignment on merged OutputHits.
// - hits: Stage 2 results (modified in-place with alignment data)
// - queries: original FASTA query sequences
// - db_path: BLAST DB path for subject sequence retrieval
// - context_is_ratio/context_ratio/context_abs: -context_extend option values
// Returns filtered hits (min_ppositive/min_npositive applied).
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
