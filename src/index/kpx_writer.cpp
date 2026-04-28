#include "index/kpx_writer.hpp"
#include "index/kpx_format.hpp"
#include "index/pfd_codec.hpp"
#include "core/config.hpp"

#include <cstdio>
#include <cstring>

namespace ikafssn {

KpxWriter::KpxWriter(int k, uint32_t freq_threshold_part)
    : k_(k),
      table_size_(ikafssn::table_size(k)),
      freq_threshold_part_(freq_threshold_part) {
    pos_offsets_.resize(table_size_, 0);
}

void KpxWriter::add_posting_list(uint32_t kmer_value,
                                  const std::vector<PostingEntry>& entries) {
    pos_offsets_[kmer_value] = posting_data_.size();
    total_postings_ += entries.size();

    if (entries.empty()) return;

    // v6: hand both seq_ids and absolute positions to the codec; entries
    // are already grouped by seq_id (caller contract).
    std::vector<uint32_t> sids(entries.size());
    std::vector<uint32_t> abs_positions(entries.size());
    for (size_t i = 0; i < entries.size(); i++) {
        sids[i]          = entries[i].seq_id;
        abs_positions[i] = entries[i].pos;
    }
    pfd::encode_posting_kpx(sids.data(), abs_positions.data(),
                            static_cast<uint32_t>(abs_positions.size()),
                            freq_threshold_part_,
                            posting_data_);
}

bool KpxWriter::write(const std::string& path) const {
    bool use_offset32 = (posting_data_.size() <= UINT32_MAX);

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
    hdr.offset_type = use_offset32 ? 0 : 1;

    // Reserved codec-extension area (codec selection follows
    // format_version since Phase 5g-1).
    hdr.codec_id      = 0;
    hdr.codec_version = 0;
    hdr.block_size    = 0;
    hdr.tail_codec    = 0;

    std::fwrite(&hdr, sizeof(hdr), 1, fp);

    // Write pos_offsets table
    if (use_offset32) {
        std::vector<uint32_t> offsets32(table_size_);
        for (uint32_t i = 0; i < table_size_; i++) {
            offsets32[i] = static_cast<uint32_t>(pos_offsets_[i]);
        }
        std::fwrite(offsets32.data(), sizeof(uint32_t), table_size_, fp);
    } else {
        std::fwrite(pos_offsets_.data(), sizeof(uint64_t), table_size_, fp);
    }

    // Write position posting data
    if (!posting_data_.empty()) {
        std::fwrite(posting_data_.data(), 1, posting_data_.size(), fp);
    }

    std::fclose(fp);
    return true;
}

} // namespace ikafssn
