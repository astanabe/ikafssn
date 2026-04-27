#include "index/kpx_writer.hpp"
#include "index/kpx_format.hpp"
#include "index/pfd_codec.hpp"
#include "core/config.hpp"

#include <cstdio>
#include <cstring>

namespace ikafssn {

KpxWriter::KpxWriter(int k)
    : k_(k), table_size_(ikafssn::table_size(k)) {
    pos_offsets_.resize(table_size_, 0);
}

void KpxWriter::add_posting_list(uint32_t kmer_value,
                                  const std::vector<PostingEntry>& entries) {
    pos_offsets_[kmer_value] = posting_data_.size();
    total_postings_ += entries.size();

    if (entries.empty()) return;

    // v4: encode absolute positions via FastPFor — sequence-boundary delta
    // reset is no longer needed because absolute positions naturally fit
    // FastPFor's per-block bit-width adaptation.
    std::vector<uint32_t> abs_positions(entries.size());
    for (size_t i = 0; i < entries.size(); i++) {
        abs_positions[i] = entries[i].pos;
    }
    pfd::encode_posting_kpx(abs_positions.data(),
                            static_cast<uint32_t>(abs_positions.size()),
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

    // v4 codec fields
    hdr.codec_id      = KPX_CODEC_PFOR_FOR;
    hdr.codec_version = 1;
    hdr.block_size    = pfd::kPfdBlockSize;
    hdr.tail_codec    = KPX_TAIL_VBYTE;

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
