#include "index/kix_writer.hpp"
#include "index/kix_format.hpp"
#include "core/config.hpp"
#include "core/varint.hpp"

#include <cstdio>
#include <cstring>

namespace ikafssn {

KixWriter::KixWriter(int k, uint8_t kmer_type)
    : k_(k), kmer_type_(kmer_type), table_size_(ikafssn::table_size(k)) {
    offsets_.resize(table_size_, 0);
    counts_.resize(table_size_, 0);
}

void KixWriter::set_volume_info(uint16_t volume_index, uint16_t total_volumes) {
    volume_index_ = volume_index;
    total_volumes_ = total_volumes;
}

void KixWriter::set_db(const std::string& name) {
    db_ = name;
}

void KixWriter::set_num_sequences(uint32_t n) {
    num_sequences_ = n;
}

void KixWriter::set_flags(uint32_t flags) {
    flags_ = flags;
}

void KixWriter::add_posting_list(uint64_t kmer_value, const std::vector<uint32_t>& seq_ids) {
    offsets_[kmer_value] = posting_data_.size();
    counts_[kmer_value] = static_cast<uint32_t>(seq_ids.size());
    total_postings_ += seq_ids.size();

    if (seq_ids.empty()) return;

    uint8_t buf[5];
    // First ID: raw varint
    size_t n = varint_encode(seq_ids[0], buf);
    posting_data_.insert(posting_data_.end(), buf, buf + n);

    // Subsequent: delta varint
    for (size_t i = 1; i < seq_ids.size(); i++) {
        uint32_t delta = seq_ids[i] - seq_ids[i - 1];
        n = varint_encode(delta, buf);
        posting_data_.insert(posting_data_.end(), buf, buf + n);
    }
}

bool KixWriter::write(const std::string& path) {
    FILE* fp = std::fopen(path.c_str(), "wb");
    if (!fp) {
        std::fprintf(stderr, "KixWriter: cannot open '%s' for writing\n", path.c_str());
        return false;
    }

    // Write header
    KixHeader hdr{};
    std::memcpy(hdr.magic, KIX_MAGIC, 4);
    hdr.format_version = KIX_FORMAT_VERSION;
    hdr.k = static_cast<uint8_t>(k_);
    hdr.kmer_type = kmer_type_;
    hdr.num_sequences = num_sequences_;
    hdr.total_postings = total_postings_;
    hdr.flags = flags_;
    hdr.volume_index = volume_index_;
    hdr.total_volumes = total_volumes_;

    size_t name_len = std::min(db_.size(), size_t(32));
    hdr.db_len = static_cast<uint16_t>(name_len);
    std::memcpy(hdr.db, db_.c_str(), name_len);

    std::fwrite(&hdr, sizeof(hdr), 1, fp);

    // Write offsets table
    std::fwrite(offsets_.data(), sizeof(uint64_t), table_size_, fp);

    // Write counts table
    std::fwrite(counts_.data(), sizeof(uint32_t), table_size_, fp);

    // Write posting data
    if (!posting_data_.empty()) {
        std::fwrite(posting_data_.data(), 1, posting_data_.size(), fp);
    }

    std::fclose(fp);
    return true;
}

} // namespace ikafssn
