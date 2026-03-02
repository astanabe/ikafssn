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

    if (std::memcmp(header_->magic, KPX_MAGIC, 4) != 0) {
        std::fprintf(stderr, "KpxReader: invalid magic\n");
        close();
        return false;
    }

    if (header_->format_version != KPX_FORMAT_VERSION) {
        std::fprintf(stderr, "KpxReader: unsupported format version %u\n",
                     header_->format_version);
        close();
        return false;
    }

    table_size_ = ikafssn::table_size(header_->k);

    const uint8_t* ptr = mmap_.data() + sizeof(KpxHeader);

    pos_offsets_ = reinterpret_cast<const uint64_t*>(ptr);
    ptr += sizeof(uint64_t) * table_size_;

    posting_data_ = ptr;
    posting_data_size_ = mmap_.size() - (ptr - mmap_.data());

    return true;
}

void KpxReader::close() {
    mmap_.close();
    header_ = nullptr;
    pos_offsets_ = nullptr;
    posting_data_ = nullptr;
    posting_data_size_ = 0;
    table_size_ = 0;
}

size_t KpxReader::willneed_size() const {
    if (!mmap_.is_open()) return 0;
    return sizeof(KpxHeader) + sizeof(uint64_t) * table_size_;
}

void KpxReader::apply_madvise(bool willneed) {
    if (!mmap_.is_open()) return;
    size_t dict_size = willneed_size();
    if (willneed) {
        mmap_.advise(0, dict_size, MADV_WILLNEED);
        mmap_.advise(dict_size, posting_data_size_, MADV_RANDOM);
    } else {
        mmap_.advise(MADV_RANDOM);
    }
    mmap_.advise(0, dict_size, MADV_HUGEPAGE);
}

} // namespace ikafssn
