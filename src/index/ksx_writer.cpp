#include "index/ksx_writer.hpp"
#include "index/ksx_format.hpp"
#include "core/config.hpp"

#include <cstdio>
#include <cstring>

namespace ikafssn {

void KsxWriter::add_sequence(uint32_t seq_length, const std::string& accession) {
    seq_lengths_.push_back(seq_length);
    accessions_.push_back(accession);
}

bool KsxWriter::write(const std::string& path) const {
    FILE* fp = std::fopen(path.c_str(), "wb");
    if (!fp) {
        std::fprintf(stderr, "KsxWriter: cannot open '%s' for writing\n", path.c_str());
        return false;
    }

    uint32_t num_seq = num_sequences();

    // Build accession offset table.  The .ksx wire format stores offsets as
    // u32, so the concatenated accession byte length must fit in 32 bits.
    // We accumulate in u64 first and abort on overflow rather than letting
    // the offsets wrap silently.
    std::vector<uint32_t> acc_offsets(num_seq + 1);
    uint64_t offset = 0;
    for (uint32_t i = 0; i < num_seq; i++) {
        if (offset > UINT32_MAX) {
            std::fprintf(stderr,
                "KsxWriter: total accession bytes exceed uint32_t at OID %u "
                "(offset=%lu).  The .ksx layout caps the accession string "
                "table at 4 GiB; split the BLAST DB into smaller volumes.\n",
                i, static_cast<unsigned long>(offset));
            std::fclose(fp);
            std::remove(path.c_str());
            return false;
        }
        acc_offsets[i] = static_cast<uint32_t>(offset);
        offset += accessions_[i].size();
    }
    if (offset > UINT32_MAX) {
        std::fprintf(stderr,
            "KsxWriter: total accession bytes exceed uint32_t (sentinel=%lu).\n",
            static_cast<unsigned long>(offset));
        std::fclose(fp);
        std::remove(path.c_str());
        return false;
    }
    acc_offsets[num_seq] = static_cast<uint32_t>(offset);

    // Write header
    KsxHeader hdr{};
    std::memcpy(hdr.magic, KSX_MAGIC, 4);
    hdr.format_version = KSX_FORMAT_VERSION;
    hdr.num_sequences = num_seq;
    std::fwrite(&hdr, sizeof(hdr), 1, fp);

    // Write seq_lengths
    std::fwrite(seq_lengths_.data(), sizeof(uint32_t), num_seq, fp);

    // Write accession_offsets (num_seq + 1 entries)
    std::fwrite(acc_offsets.data(), sizeof(uint32_t), num_seq + 1, fp);

    // Write accession strings (concatenated, no NUL terminators)
    for (uint32_t i = 0; i < num_seq; i++) {
        if (!accessions_[i].empty()) {
            std::fwrite(accessions_[i].data(), 1, accessions_[i].size(), fp);
        }
    }

    std::fclose(fp);
    return true;
}

} // namespace ikafssn
