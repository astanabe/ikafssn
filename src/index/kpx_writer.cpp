#include "index/kpx_writer.hpp"
#include "index/kpx_format.hpp"
#include "core/config.hpp"
#include "core/varint.hpp"

#include <cstdio>
#include <cstring>

namespace ikafssn {

KpxWriter::KpxWriter(int k)
    : k_(k), table_size_(ikafssn::table_size(k)) {
    pos_offsets_.resize(table_size_, 0);
}

void KpxWriter::add_posting_list(uint64_t kmer_value,
                                  const std::vector<PostingEntry>& entries) {
    pos_offsets_[kmer_value] = posting_data_.size();
    total_postings_ += entries.size();

    if (entries.empty()) return;

    uint8_t buf[5];

    // First entry: always raw pos
    size_t n = varint_encode(entries[0].pos, buf);
    posting_data_.insert(posting_data_.end(), buf, buf + n);

    // Subsequent entries
    for (size_t i = 1; i < entries.size(); i++) {
        bool new_seq = (entries[i].seq_id != entries[i - 1].seq_id);
        if (new_seq) {
            // Sequence boundary: write raw pos (delta reset)
            n = varint_encode(entries[i].pos, buf);
        } else {
            // Same sequence: write pos delta
            uint32_t delta = entries[i].pos - entries[i - 1].pos;
            n = varint_encode(delta, buf);
        }
        posting_data_.insert(posting_data_.end(), buf, buf + n);
    }
}

bool KpxWriter::write(const std::string& path) const {
    FILE* fp = std::fopen(path.c_str(), "wb");
    if (!fp) {
        std::fprintf(stderr, "KpxWriter: cannot open '%s' for writing\n", path.c_str());
        return false;
    }

    // Write header
    KpxHeader hdr{};
    std::memcpy(hdr.magic, KPX_MAGIC, 4);
    hdr.format_version = KPX_FORMAT_VERSION;
    hdr.k = static_cast<uint8_t>(k_);
    hdr.total_postings = total_postings_;

    std::fwrite(&hdr, sizeof(hdr), 1, fp);

    // Write pos_offsets table
    std::fwrite(pos_offsets_.data(), sizeof(uint64_t), table_size_, fp);

    // Write position posting data
    if (!posting_data_.empty()) {
        std::fwrite(posting_data_.data(), 1, posting_data_.size(), fp);
    }

    std::fclose(fp);
    return true;
}

} // namespace ikafssn
