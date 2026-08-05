#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "io/result_writer.hpp"
#include "io/fasta_reader.hpp"
#include "util/logger.hpp"

namespace ikafssn {

class BlastDbReader;

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
    // oversize group is processed as its own solo batch.
    uint64_t posting_budget = 0;
    // Candidate-count limits applied by the caller after run_stage3() and the
    // Stage 3 dedup, scored by alnscore.  They refine the result set but do not
    // reduce the alignment work (every chain is aligned first); cap the chains
    // entering Stage 3 via the Stage 2 limits to bound the alignment cost.
    //   N: max hits per (qseqid, sseqid[, sstrand]) group (0 = unlimited).
    //   X: selection mode 1..4 (tie / strand semantics, as in Stage 1/2).
    //   L: max hits per query, tie-inclusive (0 = unlimited).
    uint32_t max_nhit_per_subject = 1;
    uint8_t  max_nhit_per_subject_mode = 3;
    uint32_t max_nhit_in_total = 0;
};

// Run Stage 3 alignment on merged OutputHits.
// - hits: Stage 2 results (modified in-place with alignment data).  Each hit's
//   query_idx must index `queries`; hits whose query_idx or volume is out of
//   range are dropped with a warning.
// - queries: original FASTA query sequences
// - readers: one open BLAST DB reader per volume, indexed by hit.volume.
//   Owned by the caller, which lets a long-lived process keep the volumes
//   open across calls instead of re-opening them every time.
// - context_is_ratio/context_ratio/context_abs: -context_extend option values
// Returns filtered hits (min_ppositive/min_npositive applied).
std::vector<OutputHit> run_stage3(
    std::vector<OutputHit>& hits,
    const std::vector<FastaRecord>& queries,
    const std::vector<BlastDbReader>& readers,
    const Stage3Config& config,
    bool context_is_ratio,
    double context_ratio,
    uint32_t context_abs,
    const Logger& logger);

// Same, opening every volume of `db_path` for the duration of the call.
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
