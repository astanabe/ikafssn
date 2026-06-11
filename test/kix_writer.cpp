#include "kix_writer.hpp"
#include "index/kix_format.hpp"
#include "index/dictionary_io.hpp"
#include "index/pfd_codec.hpp"
#include "core/config.hpp"

#include <cstdio>
#include <cstring>

namespace ikafssn {

KixWriter::KixWriter(int k, uint8_t kmer_type)
    : k_(k), kmer_type_(kmer_type), table_size_(ikafssn::table_size(k)) {
    offsets_.resize(table_size_ + 1, 0);
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

void KixWriter::add_posting_list(uint32_t kmer_value, const std::vector<uint32_t>& seq_ids) {
    offsets_[kmer_value] = posting_file_.size();
    total_distinct_postings_ += seq_ids.size();

    if (seq_ids.empty()) return;

    // Encode the distinct seq_id delta-stream (first absolute, then
    // strictly-positive differences) via FastPFor.  Caller contract: seq_ids
    // are sorted and contain no duplicates.
    std::vector<uint32_t> deltas(seq_ids.size());
    deltas[0] = seq_ids[0];
    for (size_t i = 1; i < seq_ids.size(); i++) {
        deltas[i] = seq_ids[i] - seq_ids[i - 1];
    }
    pfd::encode_posting_kix(deltas.data(),
                            static_cast<uint32_t>(deltas.size()),
                            posting_file_);
}

bool KixWriter::write(const std::string& path) {
    offsets_[table_size_] = posting_file_.size();  // sentinel

    FILE* fp = std::fopen(path.c_str(), "wb");
    if (!fp) {
        std::fprintf(stderr, "KixWriter: cannot open '%s' for writing\n", path.c_str());
        return false;
    }

    // KIX_FLAG_OFFSET32 is reserved and forced to 0.
    KixHeader hdr{};
    std::memcpy(hdr.magic, KIX_MAGIC, sizeof(KIX_MAGIC));
    hdr.format_version = KIX_FORMAT_VERSION;
    hdr.k = static_cast<uint8_t>(k_);
    hdr.kmer_type = kmer_type_;
    hdr.num_sequences = num_sequences_;
    hdr.total_distinct_postings = total_distinct_postings_;
    hdr.flags = flags_ & ~KIX_FLAG_OFFSET32;
    hdr.volume_index = volume_index_;
    hdr.total_volumes = total_volumes_;

    size_t name_len = std::min(db_.size(), size_t(32));
    hdr.db_len = static_cast<uint16_t>(name_len);
    std::memcpy(hdr.db, db_.c_str(), name_len);

    // Reserved codec-extension area (codec selection follows
    // format_version; these fields are reserved).
    hdr.codec_id              = 0;
    hdr.codec_version         = 0;
    hdr.block_size            = 0;
    hdr.tail_codec            = 0;
    hdr.exception_codec_flags = 0;

    std::fwrite(&hdr, sizeof(hdr), 1, fp);

    if (!write_kix_dictionary_ef(fp, offsets_.data(), table_size_,
                                 posting_file_.size())) {
        std::fprintf(stderr, "KixWriter: failed to write EF dictionary to '%s'\n",
                     path.c_str());
        std::fclose(fp);
        return false;
    }

    if (!posting_file_.empty()) {
        std::fwrite(posting_file_.data(), 1, posting_file_.size(), fp);
    }

    std::fclose(fp);
    return true;
}

} // namespace ikafssn
