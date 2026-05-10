#include "ikafssnretrieve/local_retriever.hpp"
#include "io/blastdb_reader.hpp"
#include "io/accession_utils.hpp"

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

// Pick the first accession token from an sseqid that may carry a
// multi-defline '\x01'-joined string.  writes the parent
// OID's accession into hit.sseqid verbatim, so when the BLAST DB volume
// was built with `makeblastdb -parse_seqids` the field can still be a
// joined form.  We use the first token for the FASTA defline so the
// kafsss-style `parent_acc:start-end` header stays well-formed.
static std::string first_accession_token(const std::string& sseqid) {
    auto tokens = split_accessions(sseqid);
    if (tokens.empty()) return sseqid;
    return std::string(tokens.front());
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
            if (acc.empty()) continue;
            // Multi-defline OIDs register every individual accession so
            // a hit row's sseqid (which may be any one of the joined
            // accessions or the full '\x01'-joined form) resolves back
            // to the OID.  Also register the original joined form so
            // callers that pass it through verbatim still match.
            for (auto token : split_accessions(acc)) {
                acc_map[std::string(token)] = {vi, oid};
            }
            acc_map[acc] = {vi, oid};
        }
    }

    uint32_t retrieved = 0;
    for (const auto& hit : hits) {
        // hit.sseqid is the parent OID accession (Stage 2 dedup
        // collapses every fragment of one parent into a single row),
        // hit.sstart / hit.send are parent-OID-relative coordinates,
        // and BlastDbReader OIDs index parent sequences directly, so
        // hit.sseqid -> (reader, parent OID) resolves through the same
        // accession lookup with no fragment-table indirection.
        auto it = acc_map.find(hit.sseqid);
        if (it == acc_map.end()) {
            std::fprintf(stderr, "retrieve_local: accession '%s' not found in DB\n",
                         hit.sseqid.c_str());
            continue;
        }

        size_t reader_idx = it->second.first;
        uint32_t oid = it->second.second;
        uint32_t seq_len = readers[reader_idx].seq_length(oid);
        if (seq_len == 0) {
            std::fprintf(stderr, "retrieve_local: zero-length parent OID %u for '%s'\n",
                         oid, hit.sseqid.c_str());
            continue;
        }

        // Compute extraction range with context.  Coordinates are
        // 0-based parent-relative inclusive throughout this file, matching
        // the convention Stage 2/3 use for hit.sstart / hit.send.
        uint32_t ext_start = hit.sstart;
        uint32_t ext_end = hit.send;
        if (opts.context > 0) {
            ext_start = (ext_start >= opts.context) ? ext_start - opts.context : 0;
            ext_end = std::min(ext_end + opts.context, seq_len - 1);
        }
        // Clamp to actual sequence length (defensive against bad input rows).
        if (ext_end >= seq_len) {
            ext_end = seq_len - 1;
        }
        if (ext_start > ext_end) {
            std::fprintf(stderr, "retrieve_local: invalid range [%u, %u] for '%s'\n",
                         ext_start, ext_end, hit.sseqid.c_str());
            continue;
        }

        // pull only the requested window from the parent
        // OID via partial decode instead of materialising the whole
        // sequence.  On chromosome-scale parents the fragment index
        // shrinks the search window to ~min_length_split, and a full
        // get_sequence() decode would defeat that benefit.
        std::string subseq = readers[reader_idx].get_subsequence(oid, ext_start, ext_end);
        if (subseq.empty()) {
            std::fprintf(stderr, "retrieve_local: failed to get subsequence for OID %u "
                         "[%u, %u]\n", oid, ext_start, ext_end);
            continue;
        }

        // Apply reverse complement for minus strand
        if (hit.sstrand == '-') {
            reverse_complement(subseq);
        }

        // FASTA defline (kafsss style): the parent accession
        // and the (1-based inclusive) extracted range are joined into one
        // canonical sseqid token, mirroring kafssstore's
        // `accession:start:end` fragment naming and BLAST's
        // `accession:start-end` subseq notation.  Downstream tools that
        // index FASTA by the leading defline token therefore see one
        // unique sseqid per retrieved hit even when the same parent
        // accession appears in many rows.
        const std::string parent_acc = first_accession_token(hit.sseqid);
        const uint32_t one_based_start = ext_start + 1;
        const uint32_t one_based_end   = ext_end + 1;
        out << '>' << parent_acc << ':' << one_based_start << '-' << one_based_end
            << " query=" << hit.qseqid
            << " strand=" << hit.sstrand
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
