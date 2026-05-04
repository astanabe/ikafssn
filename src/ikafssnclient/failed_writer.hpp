#pragma once

#include <string>
#include <vector>

#include "io/result_writer.hpp"

namespace ikafssn {

// Synthesize per-query OutputHit rows for a failed async job.  One row
// per defline is appended to `out`, with skip_reason=kFailHttpJob and
// skip_detail=`reason`.  Downstream result_writer / sam_writer pick up
// the kFailHttpJob branch and emit *FAILED:<reason> sentinels.
//
// `deflines` is the cached defline list for the job (loaded from
// `<job_id>.deflines.zst` so we do not need the original FASTA on disk).
// Each entry should be the parsed FASTA defline up to the first
// whitespace — i.e. the qseqid that the original submit POST used.
void synth_failed_hits(const std::vector<std::string>& deflines,
                       const std::string& reason,
                       std::vector<OutputHit>& out);

} // namespace ikafssn
