#include "index/ksx_writer.hpp"
#include "index/ksx_format.hpp"
#include "core/config.hpp"

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
            std::fclose(fp);
            std::remove(path.c_str());
            return false;
        }
        acc_offsets[i] = static_cast<uint32_t>(offset);
        offset += parent_accessions_[i].size();
    }
    if (offset > UINT32_MAX) {
        std::fprintf(stderr,
            "KsxWriter: total accession bytes exceed uint32_t (sentinel=%lu).\n",
            static_cast<unsigned long>(offset));
        std::fclose(fp);
        std::remove(path.c_str());
        return false;
    }
    acc_offsets[num_parents] = static_cast<uint32_t>(offset);

    // Header
    KsxHeader hdr{};
    std::memcpy(hdr.magic, KSX_MAGIC, sizeof(KSX_MAGIC));
    hdr.format_version   = KSX_FORMAT_VERSION;
    hdr.num_sequences    = num_seq;
    hdr.num_parents      = num_parents;
    std::fwrite(&hdr, sizeof(hdr), 1, fp);

    // Parent table
    std::fwrite(parent_lengths_.data(),    sizeof(uint32_t), num_parents,     fp);
    std::fwrite(parent_blast_oids_.data(), sizeof(uint32_t), num_parents,     fp);
    std::fwrite(acc_offsets.data(),        sizeof(uint32_t), num_parents + 1, fp);
    for (uint32_t i = 0; i < num_parents; i++) {
        if (!parent_accessions_[i].empty()) {
            std::fwrite(parent_accessions_[i].data(), 1,
                        parent_accessions_[i].size(), fp);
        }
    }

    // Fragment table
    std::fwrite(fragment_parent_idx_.data(), sizeof(uint32_t), num_seq, fp);
    std::fwrite(fragment_start_.data(),      sizeof(uint32_t), num_seq, fp);
    std::fwrite(fragment_end_.data(),        sizeof(uint32_t), num_seq, fp);

    std::fclose(fp);
    return true;
}

} // namespace ikafssn
