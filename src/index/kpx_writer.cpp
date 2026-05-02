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
    pos_offsets_[kmer_value] = posting_file_.size();
    total_position_count_ += entries.size();

    if (entries.empty()) {
        // Even an empty posting list needs a valid header so the decoder can be
        // invoked with a nonzero byte length.  Emit the all-zero blob.
        pfd::encode_posting_kpx(nullptr, nullptr, 0,
                                nullptr, 0,
                                freq_threshold_part_,
                                posting_file_);
        return;
    }

    // v8: hand the encoder pre-deduplicated distinct_sid + occ_count
    // arrays alongside the sorted abs_pos array.  occ_count is u32 so a
    // single (k-mer, seq_id) cluster may exceed 255 occurrences (this
    // routinely happens in nt-class BLAST DBs with large genomic
    // contigs).  Entries are sorted by (seq_id, pos) — see header
    // contract.
    std::vector<uint32_t> distinct_sid;
    std::vector<uint32_t> occ_count;
    std::vector<uint32_t> abs_positions;
    distinct_sid.reserve(entries.size());
    occ_count.reserve(entries.size());
    abs_positions.reserve(entries.size());

    uint32_t i = 0;
    while (i < entries.size()) {
        uint32_t j = i + 1;
        while (j < entries.size() && entries[j].seq_id == entries[i].seq_id) j++;
        distinct_sid.push_back(entries[i].seq_id);
        occ_count.push_back(j - i);
        for (uint32_t e = i; e < j; e++) abs_positions.push_back(entries[e].pos);
        i = j;
    }

    pfd::encode_posting_kpx(distinct_sid.data(), occ_count.data(),
                            static_cast<uint32_t>(distinct_sid.size()),
                            abs_positions.data(),
                            static_cast<uint32_t>(abs_positions.size()),
                            freq_threshold_part_,
                            posting_file_);
}

bool KpxWriter::write(const std::string& path) const {
    bool use_offset32 = (posting_file_.size() <= UINT32_MAX);

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
    hdr.total_position_count = total_position_count_;
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

    // Write position posting file
    if (!posting_file_.empty()) {
        std::fwrite(posting_file_.data(), 1, posting_file_.size(), fp);
    }

    std::fclose(fp);
    return true;
}

} // namespace ikafssn
