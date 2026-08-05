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
        std::string acc;
        for (uint32_t oid = 0; oid < nseqs; oid++) {
            if (!readers[vi].get_accession(oid, acc)) {
                std::fprintf(stderr,
                    "retrieve_local: cannot read the accession of OID %u in "
                    "volume '%s'\n", oid, vol_paths[vi].c_str());
                return 0;
            }
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
        // hit.sstart / hit.send are parent-relative coordinates,
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

        // Per-hit context: ratio mode derives it from the hit's query length.
        // Coordinates are 0-based parent-relative half-open (the convention
        // Stage 2/3 use for hit.sstart / hit.send).
        uint32_t ctx = opts.is_ratio
            ? static_cast<uint32_t>(hit.qlen * opts.ratio)
            : opts.context;

        // Pull only the requested range from the parent OID via partial decode;
        // a full get_sequence() decode would defeat the fragment index's
        // benefit on chromosome-scale parents.
        ContextSubseq cs = extract_context_subseq(
            readers[reader_idx], oid, hit.sstart, hit.send, ctx);
        if (cs.seq_len == 0) {
            std::fprintf(stderr, "retrieve_local: zero-length parent OID %u for '%s'\n",
                         oid, hit.sseqid.c_str());
            continue;
        }
        if (cs.seq.empty()) {
            std::fprintf(stderr, "retrieve_local: failed to get subsequence for OID %u "
                         "[%u, %u)\n", oid, cs.ext_start, cs.ext_end);
            continue;
        }
        std::string subseq = std::move(cs.seq);
        uint32_t ext_start = cs.ext_start;
        uint32_t ext_end = cs.ext_end;

        // Apply reverse complement for minus strand
        if (hit.sstrand == '-') {
            reverse_complement(subseq);
        }

        // kafsss-style defline: parent_acc:start-end (1-based inclusive), so each
        // retrieved hit gets a unique leading sseqid token even when one parent
        // accession appears in many rows.
        const std::string parent_acc = first_accession_token(hit.sseqid);
        // A 0-based exclusive end and a 1-based inclusive end are the same.
        const uint32_t one_based_start = ext_start + 1;
        const uint32_t one_based_end   = ext_end;
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
