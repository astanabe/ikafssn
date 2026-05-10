#include "index/ksx_reader.hpp"
#include "index/ksx_format.hpp"
#include "core/config.hpp"

#include <cstring>
#include <cstdio>

namespace ikafssn {

bool KsxReader::open(const std::string& path) {
    close();

    if (!mmap_.open(path))
        return false;

    if (mmap_.size() < sizeof(KsxHeader)) {
        std::fprintf(stderr, "KsxReader: file too small for header\n");
        close();
        return false;
    }

    const auto* hdr = reinterpret_cast<const KsxHeader*>(mmap_.data());

    if (std::memcmp(hdr->magic, KSX_MAGIC, sizeof(KSX_MAGIC)) != 0) {
        std::fprintf(stderr, "KsxReader: invalid magic\n");
        close();
        return false;
    }

    if (hdr->format_version != KSX_FORMAT_VERSION) {
        std::fprintf(stderr,
            "KsxReader: index format version mismatch (got %u, expected %u). "
            "Please rebuild with the current ikafssnindex.\n",
            hdr->format_version, KSX_FORMAT_VERSION);
        close();
        return false;
    }

    num_sequences_    = hdr->num_sequences;
    num_parents_      = hdr->num_parents;

    const uint8_t* ptr = mmap_.data() + sizeof(KsxHeader);

    parent_lengths_ = reinterpret_cast<const uint32_t*>(ptr);
    ptr += sizeof(uint32_t) * num_parents_;

    parent_blast_oids_ = reinterpret_cast<const uint32_t*>(ptr);
    ptr += sizeof(uint32_t) * num_parents_;

    parent_acc_offsets_ = reinterpret_cast<const uint32_t*>(ptr);
    ptr += sizeof(uint32_t) * (num_parents_ + 1);

    parent_acc_strings_ = reinterpret_cast<const char*>(ptr);
    if (num_parents_ > 0) {
        ptr += parent_acc_offsets_[num_parents_];
    }

    fragment_parent_idx_ = reinterpret_cast<const uint32_t*>(ptr);
    ptr += sizeof(uint32_t) * num_sequences_;

    fragment_start_ = reinterpret_cast<const uint32_t*>(ptr);
    ptr += sizeof(uint32_t) * num_sequences_;

    fragment_end_ = reinterpret_cast<const uint32_t*>(ptr);
    ptr += sizeof(uint32_t) * num_sequences_;

    if (static_cast<size_t>(ptr - mmap_.data()) > mmap_.size()) {
        std::fprintf(stderr, "KsxReader: file truncated\n");
        close();
        return false;
    }

    return true;
}

void KsxReader::close() {
    mmap_.close();
    num_sequences_   = 0;
    num_parents_     = 0;
    parent_lengths_      = nullptr;
    parent_blast_oids_   = nullptr;
    parent_acc_offsets_  = nullptr;
    parent_acc_strings_  = nullptr;
    fragment_parent_idx_ = nullptr;
    fragment_start_      = nullptr;
    fragment_end_        = nullptr;
}

uint32_t KsxReader::seq_length(uint32_t seq_id) const {
    return fragment_end_[seq_id] - fragment_start_[seq_id] + 1u;
}

std::string_view KsxReader::accession(uint32_t seq_id) const {
    return parent_accession(fragment_parent_idx_[seq_id]);
}

std::string_view KsxReader::parent_accession(uint32_t parent_idx) const {
    uint32_t start = parent_acc_offsets_[parent_idx];
    uint32_t len   = parent_acc_offsets_[parent_idx + 1] - start;
    return std::string_view(parent_acc_strings_ + start, len);
}

size_t KsxReader::willneed_size() const {
    if (!mmap_.is_open()) return 0;
    return mmap_.size();
}

void KsxReader::apply_madvise(bool willneed) {
    mmap_.advise_all(willneed);
}

} // namespace ikafssn
