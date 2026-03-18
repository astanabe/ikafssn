#include "io/primer_query.hpp"
#include "core/spaced_seed.hpp"
#include "core/kmer_encoding.hpp"

#include <set>
#include <cstdio>

namespace ikafssn {

// Count k-mer positions in a sequence using KmerScanner.
// For spaced seeds (t > 0), uses scan_spaced; otherwise uses scan.
static uint32_t count_kmer_positions(const std::string& seq, int k,
                                      uint8_t t,
                                      const std::vector<uint32_t>* masks) {
    std::set<uint32_t> positions;

    if (k <= 8) {
        KmerScanner<uint16_t> scanner(k);
        if (t > 0 && masks && !masks->empty()) {
            scanner.scan_spaced(seq.c_str(), seq.size(), *masks, t,
                [&positions](uint32_t pos, uint16_t) { positions.insert(pos); });
        } else {
            scanner.scan(seq.c_str(), seq.size(),
                [&positions](uint32_t pos, uint16_t) { positions.insert(pos); });
        }
    } else {
        KmerScanner<uint32_t> scanner(k);
        if (t > 0 && masks && !masks->empty()) {
            scanner.scan_spaced(seq.c_str(), seq.size(), *masks, t,
                [&positions](uint32_t pos, uint32_t) { positions.insert(pos); });
        } else {
            scanner.scan(seq.c_str(), seq.size(),
                [&positions](uint32_t pos, uint32_t) { positions.insert(pos); });
        }
    }

    return static_cast<uint32_t>(positions.size());
}

std::string parse_primer_pairs(
    const std::vector<FastaRecord>& records,
    const PrimerConfig& config,
    std::vector<PrimerPair>& pairs) {

    pairs.clear();

    if (records.empty()) {
        return "Error: primer FASTA is empty";
    }
    if (records.size() % 2 != 0) {
        return "Error: primer FASTA must contain an even number of sequences (got " +
               std::to_string(records.size()) + ")";
    }
    if (config.insert_length < static_cast<uint32_t>(config.k)) {
        return "Error: -insert_length (" + std::to_string(config.insert_length) +
               ") must be >= k (" + std::to_string(config.k) + ")";
    }

    pairs.reserve(records.size() / 2);

    for (size_t i = 0; i < records.size(); i += 2) {
        const auto& fwd = records[i];
        const auto& rev = records[i + 1];

        if (static_cast<int>(fwd.sequence.size()) < config.k) {
            return "Error: forward primer '" + fwd.id + "' length (" +
                   std::to_string(fwd.sequence.size()) + ") < k (" +
                   std::to_string(config.k) + ")";
        }
        if (static_cast<int>(rev.sequence.size()) < config.k) {
            return "Error: reverse primer '" + rev.id + "' length (" +
                   std::to_string(rev.sequence.size()) + ") < k (" +
                   std::to_string(config.k) + ")";
        }

        PrimerPair pair;
        pair.fwd_id = fwd.id;
        pair.rev_id = rev.id;
        pair.fwd_seq = fwd.sequence;
        pair.rev_seq = rev.sequence;

        // Generate reverse complement of reverse primer
        std::string rc_rev = reverse_complement_string(rev.sequence);

        // Build query: fwd + N×insert_length + RC(rev)
        pair.query_id = fwd.id + "__+__" + rev.id;
        pair.query_seq = fwd.sequence +
                         std::string(config.insert_length, 'N') +
                         rc_rev;

        // Count k-mer positions for each primer region
        pair.fwd_kmer_positions = count_kmer_positions(
            fwd.sequence, config.k, config.t, config.masks);
        pair.rev_kmer_positions = count_kmer_positions(
            rc_rev, config.k, config.t, config.masks);

        pairs.push_back(std::move(pair));
    }

    return "";
}

} // namespace ikafssn
