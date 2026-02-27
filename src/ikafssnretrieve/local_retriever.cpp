#include "ikafssnretrieve/local_retriever.hpp"
#include "io/blastdb_reader.hpp"

#include <algorithm>
#include <cstdio>
#include <unordered_map>

namespace ikafssn {

// Reverse complement a DNA string in-place.
static void reverse_complement(std::string& seq) {
    std::reverse(seq.begin(), seq.end());
    for (auto& c : seq) {
        switch (c) {
            case 'A': c = 'T'; break;
            case 'T': c = 'A'; break;
            case 'C': c = 'G'; break;
            case 'G': c = 'C'; break;
            // N and others stay as-is
        }
    }
}

uint32_t retrieve_local(const std::vector<OutputHit>& hits,
                        const std::string& db_path,
                        const RetrieveOptions& opts,
                        std::ostream& out) {
    // Find all BLAST DB volumes
    auto vol_paths = BlastDbReader::find_volume_paths(db_path);
    if (vol_paths.empty()) {
        std::fprintf(stderr, "retrieve_local: no volumes found for DB '%s'\n",
                     db_path.c_str());
        return 0;
    }

    // Open all volumes and build accession -> (volume_index, oid) map
    std::vector<BlastDbReader> readers(vol_paths.size());
    // accession -> (reader index, oid)
    std::unordered_map<std::string, std::pair<size_t, uint32_t>> acc_map;

    for (size_t vi = 0; vi < vol_paths.size(); vi++) {
        if (!readers[vi].open(vol_paths[vi])) {
            std::fprintf(stderr, "retrieve_local: cannot open volume '%s'\n",
                         vol_paths[vi].c_str());
            return 0;
        }
        uint32_t nseqs = readers[vi].num_sequences();
        for (uint32_t oid = 0; oid < nseqs; oid++) {
            std::string acc = readers[vi].get_accession(oid);
            if (!acc.empty()) {
                acc_map[acc] = {vi, oid};
            }
        }
    }

    uint32_t retrieved = 0;
    for (const auto& hit : hits) {
        auto it = acc_map.find(hit.accession);
        if (it == acc_map.end()) {
            std::fprintf(stderr, "retrieve_local: accession '%s' not found in DB\n",
                         hit.accession.c_str());
            continue;
        }

        size_t reader_idx = it->second.first;
        uint32_t oid = it->second.second;
        uint32_t seq_len = readers[reader_idx].seq_length(oid);

        // Compute extraction range with context
        uint32_t ext_start = hit.s_start;
        uint32_t ext_end = hit.s_end;
        if (opts.context > 0) {
            ext_start = (ext_start >= opts.context) ? ext_start - opts.context : 0;
            ext_end = std::min(ext_end + opts.context, seq_len - 1);
        }

        // Get full sequence and extract substring
        std::string full_seq = readers[reader_idx].get_sequence(oid);
        if (full_seq.empty()) {
            std::fprintf(stderr, "retrieve_local: failed to get sequence for OID %u\n", oid);
            continue;
        }

        // Clamp to actual sequence length
        if (ext_end >= full_seq.size()) {
            ext_end = static_cast<uint32_t>(full_seq.size()) - 1;
        }
        if (ext_start > ext_end) {
            std::fprintf(stderr, "retrieve_local: invalid range [%u, %u] for '%s'\n",
                         ext_start, ext_end, hit.accession.c_str());
            continue;
        }

        std::string subseq = full_seq.substr(ext_start, ext_end - ext_start + 1);

        // Apply reverse complement for minus strand
        if (hit.strand == '-') {
            reverse_complement(subseq);
        }

        // Write FASTA record
        // >accession query=query_id strand=+/- range=start-end
        out << '>' << hit.accession
            << " query=" << hit.query_id
            << " strand=" << hit.strand
            << " range=" << ext_start << '-' << ext_end
            << " score=" << hit.chainscore
            << '\n';

        // Write sequence in 70-char lines
        for (size_t i = 0; i < subseq.size(); i += 70) {
            size_t len = std::min<size_t>(70, subseq.size() - i);
            out.write(subseq.data() + i, static_cast<std::streamsize>(len));
            out << '\n';
        }

        retrieved++;
    }

    return retrieved;
}

} // namespace ikafssn
