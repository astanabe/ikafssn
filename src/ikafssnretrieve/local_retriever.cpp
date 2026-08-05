#include "ikafssnretrieve/local_retriever.hpp"
#include "io/blastdb_reader.hpp"
#include "io/accession_utils.hpp"

#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include <fcntl.h>
#include <unistd.h>

namespace ikafssn {

namespace {

// Reverse complement a DNA string in-place.
void reverse_complement(std::string& seq) {
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
// multi-defline '\x01'-joined string.  The search writes the parent OID's
// accession into hit.sseqid verbatim, so when the BLAST DB volume was built
// with `makeblastdb -parse_seqids` the field can be a joined form.
std::string first_accession_token(const std::string& sseqid) {
    auto tokens = split_accessions(sseqid);
    if (tokens.empty()) return sseqid;
    return std::string(tokens.front());
}

// A scratch file the kernel reclaims however the process ends: it is
// unlinked the moment it exists, so only the descriptor keeps it alive.
class ScratchFile {
public:
    ~ScratchFile() { if (fd_ >= 0) ::close(fd_); }

    bool open(const std::string& dir) {
        std::string path = dir.empty() ? std::string(".") : dir;
        if (path.back() != '/') path.push_back('/');
        path += "ikafssnretrieve.";
        path += std::to_string(static_cast<long>(::getpid()));
        path += ".tmp";

        fd_ = ::open(path.c_str(), O_RDWR | O_CREAT | O_EXCL, 0600);
        if (fd_ < 0) {
            std::fprintf(stderr,
                "retrieve_local: cannot create the scratch file '%s': %s\n",
                path.c_str(), std::strerror(errno));
            return false;
        }
        if (::unlink(path.c_str()) != 0) {
            std::fprintf(stderr,
                "retrieve_local: cannot unlink the scratch file '%s': %s\n",
                path.c_str(), std::strerror(errno));
            ::close(fd_);
            fd_ = -1;
            return false;
        }
        return true;
    }

    // Append `data` and report where it landed.
    bool append(const std::string& data, uint64_t& offset) {
        offset = size_;
        const char* p = data.data();
        size_t left = data.size();
        while (left > 0) {
            ssize_t n = ::write(fd_, p, left);
            if (n < 0) {
                if (errno == EINTR) continue;
                std::fprintf(stderr,
                    "retrieve_local: cannot write to the scratch file: %s\n",
                    std::strerror(errno));
                return false;
            }
            p += n;
            left -= static_cast<size_t>(n);
            size_ += static_cast<uint64_t>(n);
        }
        return true;
    }

    bool read_at(uint64_t offset, size_t len, std::string& out) {
        out.resize(len);
        char* p = out.data();
        size_t left = len;
        while (left > 0) {
            ssize_t n = ::pread(fd_, p, left, static_cast<off_t>(offset));
            if (n < 0) {
                if (errno == EINTR) continue;
                std::fprintf(stderr,
                    "retrieve_local: cannot read back the scratch file: %s\n",
                    std::strerror(errno));
                return false;
            }
            if (n == 0) {
                std::fprintf(stderr,
                    "retrieve_local: the scratch file ended early\n");
                return false;
            }
            p += n;
            left -= static_cast<size_t>(n);
            offset += static_cast<uint64_t>(n);
        }
        return true;
    }

private:
    int fd_ = -1;
    uint64_t size_ = 0;
};

// One staged FASTA record, keyed by the position of its hit in the input so
// the output can be restored to input order.
struct StagedRecord {
    uint32_t hit_index;
    uint64_t offset;
    uint64_t length;
};

const char* const kAccessionLookupRequirement =
    "The BLAST DB must be searchable by accession: a BLAST v5 database "
    "(LMDB), or a v4 database built with `makeblastdb -parse_seqids`.  A "
    "database other than the one the index was built from will also fail to "
    "resolve.";

} // namespace

uint32_t retrieve_local(const std::vector<OutputHit>& hits,
                        const std::string& db_path,
                        const RetrieveOptions& opts,
                        std::ostream& out) {
    auto vol_paths = BlastDbReader::find_volume_paths(db_path);
    if (vol_paths.empty()) {
        std::fprintf(stderr, "retrieve_local: no volumes found for DB '%s'\n",
                     db_path.c_str());
        return 0;
    }

    // Group the hits by the volume the search recorded for them, so only the
    // volumes that carry a hit are ever opened.
    std::vector<std::vector<uint32_t>> hits_by_volume(vol_paths.size());
    for (uint32_t i = 0; i < hits.size(); i++) {
        if (hits[i].volume >= vol_paths.size()) {
            std::fprintf(stderr,
                "retrieve_local: hit '%s' names volume %u, but DB '%s' has "
                "only %zu volume(s)\n",
                hits[i].sseqid.c_str(), static_cast<unsigned>(hits[i].volume),
                db_path.c_str(), vol_paths.size());
            continue;
        }
        hits_by_volume[hits[i].volume].push_back(i);
    }

    ScratchFile scratch;
    if (!scratch.open(opts.scratch_dir)) return 0;

    std::vector<StagedRecord> staged;
    staged.reserve(hits.size());
    bool requirement_shown = false;

    for (size_t vi = 0; vi < vol_paths.size(); vi++) {
        if (hits_by_volume[vi].empty()) continue;

        // Exactly one volume open at a time: CSeqDB holds a descriptor per
        // mapped volume file, and the atlas owning the mappings is a
        // process-wide singleton, so descriptors return only once the last
        // reader is destroyed.
        BlastDbReader reader;
        if (!reader.open(vol_paths[vi])) {
            std::fprintf(stderr, "retrieve_local: cannot open volume '%s'\n",
                         vol_paths[vi].c_str());
            return 0;
        }

        // Resolve each hit's accession to an OID within this volume, then walk
        // the volume in OID order so the mmap is read front to back.
        struct Resolved {
            uint32_t hit_index;
            uint32_t oid;
        };
        std::vector<Resolved> resolved;
        resolved.reserve(hits_by_volume[vi].size());

        std::vector<uint32_t> oids;
        for (uint32_t hit_index : hits_by_volume[vi]) {
            const OutputHit& hit = hits[hit_index];
            const std::string acc = first_accession_token(hit.sseqid);
            if (!reader.accession_to_oids(acc, oids)) {
                std::fprintf(stderr,
                    "retrieve_local: accession lookup failed for '%s' in "
                    "volume '%s'\n", acc.c_str(), vol_paths[vi].c_str());
                if (!requirement_shown) {
                    std::fprintf(stderr, "retrieve_local: %s\n",
                                 kAccessionLookupRequirement);
                    requirement_shown = true;
                }
                continue;
            }
            if (oids.empty()) {
                std::fprintf(stderr,
                    "retrieve_local: accession '%s' not found in volume '%s'\n",
                    acc.c_str(), vol_paths[vi].c_str());
                if (!requirement_shown) {
                    std::fprintf(stderr, "retrieve_local: %s\n",
                                 kAccessionLookupRequirement);
                    requirement_shown = true;
                }
                continue;
            }
            if (oids.size() > 1) {
                // Narrow a shared accession down by the parent length the hit
                // row carries; without a unique answer the sequence to extract
                // is undetermined.
                std::vector<uint32_t> matching;
                for (uint32_t oid : oids) {
                    if (hit.slen != 0 && reader.seq_length(oid) == hit.slen)
                        matching.push_back(oid);
                }
                if (matching.size() != 1) {
                    std::fprintf(stderr,
                        "retrieve_local: accession '%s' maps to %zu OIDs in "
                        "volume '%s' and slen=%u does not single one out\n",
                        acc.c_str(), oids.size(), vol_paths[vi].c_str(),
                        hit.slen);
                    continue;
                }
                oids = matching;
            }

            resolved.push_back({hit_index, oids.front()});
        }

        std::sort(resolved.begin(), resolved.end(),
                  [](const Resolved& a, const Resolved& b) {
                      return a.oid < b.oid;
                  });

        std::string record;
        for (const Resolved& r : resolved) {
            const OutputHit& hit = hits[r.hit_index];

            // Per-hit context: ratio mode derives it from the hit's query
            // length.  Coordinates are 0-based parent-relative half-open (the
            // convention Stage 2/3 use for hit.sstart / hit.send).
            uint32_t ctx = opts.is_ratio
                ? static_cast<uint32_t>(hit.qlen * opts.ratio)
                : opts.context;

            // Pull only the requested range from the parent OID via partial
            // decode; a full get_sequence() decode would defeat the fragment
            // index's benefit on chromosome-scale parents.
            ContextSubseq cs = extract_context_subseq(
                reader, r.oid, hit.sstart, hit.send, ctx);
            if (cs.seq_len == 0) {
                std::fprintf(stderr,
                    "retrieve_local: zero-length parent OID %u for '%s'\n",
                    r.oid, hit.sseqid.c_str());
                continue;
            }
            if (cs.seq.empty()) {
                std::fprintf(stderr,
                    "retrieve_local: failed to get subsequence for OID %u "
                    "[%u, %u)\n", r.oid, cs.ext_start, cs.ext_end);
                continue;
            }
            std::string subseq = std::move(cs.seq);

            if (hit.sstrand == '-') {
                reverse_complement(subseq);
            }

            // Defline ID: `parent_acc:start-end`, 1-based inclusive, unique
            // per record even when one parent accession appears in many rows.
            // A 0-based exclusive end and a 1-based inclusive end are equal.
            // The '\x01'-joined form goes in `sseqid=` instead, since in the
            // ID it would break both the split-on-'\x01' contract in
            // io/accession_utils.hpp and any parser that reads the ID as
            // everything up to the first space.
            record.clear();
            record += '>';
            record += first_accession_token(hit.sseqid);
            record += ':';
            record += std::to_string(cs.ext_start + 1);
            record += '-';
            record += std::to_string(cs.ext_end);
            record += " sseqid=";
            record += hit.sseqid;
            record += " query=";
            record += hit.qseqid;
            record += " strand=";
            record += hit.sstrand;
            record += " score=";
            record += std::to_string(hit.chainscore);
            record += '\n';

            // Sequence in 70-char lines
            for (size_t i = 0; i < subseq.size(); i += 70) {
                size_t len = std::min<size_t>(70, subseq.size() - i);
                record.append(subseq, i, len);
                record += '\n';
            }

            uint64_t offset = 0;
            if (!scratch.append(record, offset)) return 0;
            staged.push_back({r.hit_index, offset, record.size()});
        }
    }

    // Restore input order: the volumes were walked in volume/OID order, not
    // in the order the hits arrived.
    std::sort(staged.begin(), staged.end(),
              [](const StagedRecord& a, const StagedRecord& b) {
                  return a.hit_index < b.hit_index;
              });

    uint32_t retrieved = 0;
    std::string record;
    for (const StagedRecord& s : staged) {
        if (!scratch.read_at(s.offset, static_cast<size_t>(s.length), record))
            return retrieved;
        out.write(record.data(), static_cast<std::streamsize>(record.size()));
        retrieved++;
    }

    return retrieved;
}

} // namespace ikafssn
