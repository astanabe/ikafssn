#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "io/fasta_reader.hpp"

namespace ikafssn {

struct PrimerPair {
    std::string fwd_id, rev_id;
    std::string fwd_seq, rev_seq;
    std::string query_id;       // "fwd_id__+__rev_id"
    std::string query_seq;      // fwd + N×insert_len + RC(rev)
    uint32_t fwd_kmer_positions;
    uint32_t rev_kmer_positions;
};

struct PrimerConfig {
    uint32_t insert_length;
    int k;
    uint8_t t = 0;
    const std::vector<uint32_t>* masks = nullptr;
};

// Parse FASTA records into primer pairs.
// Records must be even in number; consecutive pairs form (fwd, rev).
// Returns empty string on success, error message on failure.
std::string parse_primer_pairs(
    const std::vector<FastaRecord>& records,
    const PrimerConfig& config,
    std::vector<PrimerPair>& pairs);

} // namespace ikafssn
