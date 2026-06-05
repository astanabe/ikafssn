#include "index/kpx_reader.hpp"
#include "core/config.hpp"

#include <sys/mman.h>

#include <cstring>
#include <cstdio>

namespace ikafssn {

bool KpxReader::open(const std::string& path) {
    close();

    if (!mmap_.open(path))
        return false;

    if (mmap_.size() < sizeof(KpxHeader)) {
        std::fprintf(stderr, "KpxReader: file too small for header\n");
        close();
        return false;
    }

    header_ = reinterpret_cast<const KpxHeader*>(mmap_.data());

    if (std::memcmp(header_->magic, KPX_MAGIC, sizeof(KPX_MAGIC)) != 0) {
        std::fprintf(stderr, "KpxReader: invalid magic\n");
        close();
        return false;
    }

    if (header_->format_version != KPX_FORMAT_VERSION) {
        std::fprintf(stderr,
            "KpxReader: index format version mismatch (got %u, expected %u). "
            "Please rebuild with the current ikafssnindex.\n",
            header_->format_version, KPX_FORMAT_VERSION);
        close();
        return false;
    }

    table_size_ = ikafssn::table_size(header_->k);

    // pos_offsets dictionary is Elias-Fano; offset_type is reserved
    // and ignored at read time.
    const uint8_t* dict_ptr = mmap_.data() + sizeof(KpxHeader);
    const size_t   remaining = mmap_.size() - sizeof(KpxHeader);
    if (!pos_dict_.open(dict_ptr, remaining)) {
        std::fprintf(stderr, "KpxReader: invalid Elias-Fano dictionary\n");
        close();
        return false;
    }
    if (pos_dict_.size() != static_cast<size_t>(table_size_)) {
        std::fprintf(stderr,
            "KpxReader: pos_offsets entry count %zu does not match expected %zu\n",
            pos_dict_.size(), static_cast<size_t>(table_size_));
        close();
        return false;
    }

    const uint8_t* ptr = dict_ptr + pos_dict_.blob_bytes();
    posting_file_ = ptr;
    posting_file_size_ = mmap_.size() - (ptr - mmap_.data());

    return true;
}

void KpxReader::close() {
    mmap_.close();
    header_ = nullptr;
    pos_dict_ = ef::EFDictionary{};
    posting_file_ = nullptr;
    posting_file_size_ = 0;
    table_size_ = 0;
}

size_t KpxReader::willneed_size() const {
    if (!mmap_.is_open()) return 0;
    return sizeof(KpxHeader) + pos_dict_.willneed_size();
}

void KpxReader::apply_madvise(bool willneed) {
    mmap_.advise_dict_posting(willneed_size(), willneed);
}

size_t KpxReader::willneed_size_full() const {
    if (!mmap_.is_open()) return 0;
    return mmap_.size();
}


void KpxReader::apply_madvise_posting_random() {
    if (!mmap_.is_open()) return;
    const size_t dict_size = willneed_size();
    mmap_.advise(0, dict_size, MADV_WILLNEED);
#ifdef MADV_HUGEPAGE
    mmap_.advise(0, dict_size, MADV_HUGEPAGE);
#endif
    if (mmap_.size() > dict_size)
        mmap_.advise(dict_size, mmap_.size() - dict_size, MADV_RANDOM);
}

void KpxReader::apply_madvise_dict_only() {
    if (!mmap_.is_open()) return;
    const size_t dict = willneed_size();
    mmap_.advise(0, dict, MADV_WILLNEED);
#ifdef MADV_HUGEPAGE
    mmap_.advise(0, dict, MADV_HUGEPAGE);
#endif
    if (mmap_.size() > dict)
        mmap_.advise(dict, mmap_.size() - dict, MADV_RANDOM);
}

size_t KpxReader::dict_size() const {
    return willneed_size();
}

void KpxReader::apply_madvise_posting_ranges(
    const std::vector<std::pair<uint64_t, uint64_t>>& ranges) {
    if (!mmap_.is_open()) return;
    const size_t base = dict_size();
    for (const auto& r : ranges) {
        if (r.second == 0) continue;
        mmap_.advise(base + r.first, r.second, MADV_WILLNEED);
    }
}

} // namespace ikafssn
