#include "index/kix_reader.hpp"
#include "core/config.hpp"

#include <sys/mman.h>
#include <cstring>
#include <cstdio>

namespace ikafssn {

bool KixReader::open(const std::string& path) {
    close();

    if (!mmap_.open(path))
        return false;

    if (mmap_.size() < sizeof(KixHeader)) {
        std::fprintf(stderr, "KixReader: file too small for header\n");
        close();
        return false;
    }

    header_ = reinterpret_cast<const KixHeader*>(mmap_.data());

    if (std::memcmp(header_->magic, KIX_MAGIC, 4) != 0) {
        std::fprintf(stderr, "KixReader: invalid magic\n");
        close();
        return false;
    }

    if (header_->format_version != KIX_FORMAT_VERSION) {
        std::fprintf(stderr, "KixReader: unsupported format version %u\n",
                     header_->format_version);
        close();
        return false;
    }

    table_size_ = ikafssn::table_size(header_->k);

    const uint8_t* ptr = mmap_.data() + sizeof(KixHeader);

    offsets_ = reinterpret_cast<const uint64_t*>(ptr);
    ptr += sizeof(uint64_t) * table_size_;

    counts_ = reinterpret_cast<const uint32_t*>(ptr);
    ptr += sizeof(uint32_t) * table_size_;

    posting_data_ = ptr;
    posting_data_size_ = mmap_.size() - (ptr - mmap_.data());

    return true;
}

void KixReader::close() {
    mmap_.close();
    header_ = nullptr;
    offsets_ = nullptr;
    counts_ = nullptr;
    posting_data_ = nullptr;
    posting_data_size_ = 0;
    table_size_ = 0;
}

size_t KixReader::willneed_size() const {
    if (!mmap_.is_open()) return 0;
    return sizeof(KixHeader)
         + sizeof(uint64_t) * table_size_
         + sizeof(uint32_t) * table_size_;
}

void KixReader::apply_madvise(bool willneed) {
    if (!mmap_.is_open()) return;
    size_t dict_size = willneed_size();
    if (willneed) {
        mmap_.advise(0, dict_size, MADV_WILLNEED);
        mmap_.advise(dict_size, posting_data_size_, MADV_RANDOM);
    } else {
        mmap_.advise(MADV_RANDOM);
    }
#ifdef MADV_HUGEPAGE
    mmap_.advise(0, dict_size, MADV_HUGEPAGE);
#endif
}

} // namespace ikafssn
