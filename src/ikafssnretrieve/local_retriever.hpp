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
    // Directory for the scratch file the records are staged in.  It needs
    // room for the uncompressed output.  TMPDIR is deliberately not
    // consulted: it is a tmpfs on many systems, which would put the staged
    // records back in the memory this staging exists to save.
    std::string scratch_dir = ".";
};

// Retrieve matched subsequences from a local BLAST DB.
// Writes FASTA records to the output stream in the order the hits arrive.
// Returns the number of successfully retrieved sequences.
//
// Each hit's `volume` selects the BLAST DB volume to read and its `sseqid`
// is resolved to an OID within that volume, so the database must be
// searchable by accession (a BLAST v5 LMDB, or a v4 string ISAM built with
// `makeblastdb -parse_seqids`).  Volumes are opened one at a time and the
// records staged under `opts.scratch_dir`, which is what keeps the output
// independent of the order the volumes are walked in.
uint32_t retrieve_local(const std::vector<OutputHit>& hits,
                        const std::string& db_path,
                        const RetrieveOptions& opts,
                        std::ostream& out);

} // namespace ikafssn
