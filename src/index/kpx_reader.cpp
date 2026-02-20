#include "index/kpx_reader.hpp"
#include "core/config.hpp"

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

} // namespace ikafssn
