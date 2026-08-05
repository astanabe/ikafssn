#include "index/ksx_writer.hpp"
#include "index/ksx_format.hpp"
#include "core/config.hpp"

#include <cerrno>
#include <cstdio>
#include <cstring>

namespace ikafssn {

uint32_t KsxWriter::add_parent(uint32_t blast_oid,
                               uint32_t parent_length,
                               const std::string& accession) {
    const uint32_t idx = static_cast<uint32_t>(parent_lengths_.size());
    parent_lengths_.push_back(parent_length);
    parent_blast_oids_.push_back(blast_oid);
    parent_accessions_.push_back(accession);
    return idx;
}

uint32_t KsxWriter::add_fragment(uint32_t parent_idx,
                                 uint32_t frag_start,
                                 uint32_t frag_end) {
    const uint32_t seq_id = static_cast<uint32_t>(fragment_parent_idx_.size());
    fragment_parent_idx_.push_back(parent_idx);
    fragment_start_.push_back(frag_start);
    fragment_end_.push_back(frag_end);
    return seq_id;
}

bool KsxWriter::write(const std::string& path) const {
    FILE* fp = std::fopen(path.c_str(), "wb");
    if (!fp) {
        std::fprintf(stderr, "KsxWriter: cannot open '%s' for writing\n", path.c_str());
        return false;
    }

    const uint32_t num_parents = this->num_parents();
    const uint32_t num_seq = this->num_sequences();

    auto fail = [&]() {
        std::fclose(fp);
        std::remove(path.c_str());
        return false;
    };

    // A short `fwrite` and a failing `fclose` (the final flush) are how a
    // full disk shows up; neither is visible in the data itself.
    auto put = [&](const void* data, size_t size, size_t count) {
        if (size == 0 || count == 0) return true;
        if (std::fwrite(data, size, count, fp) != count) {
            std::fprintf(stderr,
                "KsxWriter: failed to write %lu byte(s) to '%s': %s\n",
                static_cast<unsigned long>(size * count),
                path.c_str(), std::strerror(errno));
            return false;
        }
        return true;
    };

    // Build the parent accession offset array.  The .ksx wire format stores
    // offsets as u32, so the concatenated accession byte length must fit
    // in 32 bits.  Accumulate in u64 first and abort on overflow.
    std::vector<uint32_t> acc_offsets(num_parents + 1);
    uint64_t offset = 0;
    for (uint32_t i = 0; i < num_parents; i++) {
        if (offset > UINT32_MAX) {
            std::fprintf(stderr,
                "KsxWriter: total accession bytes exceed uint32_t at parent %u "
                "(offset=%lu).  The .ksx layout caps the accession string "
                "table at 4 GiB; split the BLAST DB into smaller volumes.\n",
                i, static_cast<unsigned long>(offset));
            return fail();
        }
        acc_offsets[i] = static_cast<uint32_t>(offset);
        offset += parent_accessions_[i].size();
    }
    if (offset > UINT32_MAX) {
        std::fprintf(stderr,
            "KsxWriter: total accession bytes exceed uint32_t (sentinel=%lu).\n",
            static_cast<unsigned long>(offset));
        return fail();
    }
    acc_offsets[num_parents] = static_cast<uint32_t>(offset);

    // Header
    KsxHeader hdr{};
    std::memcpy(hdr.magic, KSX_MAGIC, sizeof(KSX_MAGIC));
    hdr.format_version   = KSX_FORMAT_VERSION;
    hdr.num_sequences    = num_seq;
    hdr.num_parents      = num_parents;
    if (!put(&hdr, sizeof(hdr), 1)) return fail();

    // Parent table
    if (!put(parent_lengths_.data(),    sizeof(uint32_t), num_parents) ||
        !put(parent_blast_oids_.data(), sizeof(uint32_t), num_parents) ||
        !put(acc_offsets.data(),        sizeof(uint32_t), num_parents + 1))
        return fail();
    for (uint32_t i = 0; i < num_parents; i++) {
        if (!put(parent_accessions_[i].data(), 1, parent_accessions_[i].size()))
            return fail();
    }

    // Fragment table
    if (!put(fragment_parent_idx_.data(), sizeof(uint32_t), num_seq) ||
        !put(fragment_start_.data(),      sizeof(uint32_t), num_seq) ||
        !put(fragment_end_.data(),        sizeof(uint32_t), num_seq))
        return fail();

    if (std::fclose(fp) != 0) {
        std::fprintf(stderr, "KsxWriter: failed to close '%s': %s\n",
                     path.c_str(), std::strerror(errno));
        std::remove(path.c_str());
        return false;
    }
    return true;
}

} // namespace ikafssn
