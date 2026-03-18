#include "index/khx_reader.hpp"
#include "index/khx_format.hpp"
#include "core/config.hpp"

#include <sys/mman.h>
#include <cstring>
#include <cstdio>

namespace ikafssn {

bool KhxReader::open(const std::string& path) {
    close();

    if (!mmap_.open(path, /*quiet=*/true))
        return false;

    if (mmap_.size() < sizeof(KhxHeader)) {
        std::fprintf(stderr, "KhxReader: file too small for header\n");
        close();
        return false;
    }

    const auto* hdr = reinterpret_cast<const KhxHeader*>(mmap_.data());

    if (std::memcmp(hdr->magic, KHX_MAGIC, 4) != 0) {
        std::fprintf(stderr, "KhxReader: invalid magic\n");
        close();
        return false;
    }

    if (hdr->format_version != KHX_FORMAT_VERSION) {
        std::fprintf(stderr, "KhxReader: unsupported format version %u\n", hdr->format_version);
        close();
        return false;
    }

    k_ = hdr->k;
    tbl_size_ = table_size(k_);
    t_ = hdr->t;
    template_type_ = hdr->template_type;
    uint32_t bitset_bytes = (tbl_size_ + 7) / 8;

    if (mmap_.size() < sizeof(KhxHeader) + bitset_bytes) {
        std::fprintf(stderr, "KhxReader: file too small for bitset data\n");
        close();
        return false;
    }

    bitset_ = mmap_.data() + sizeof(KhxHeader);

    return true;
}

void KhxReader::close() {
    mmap_.close();
    k_ = 0;
    t_ = 0;
    template_type_ = 0;
    bitset_ = nullptr;
    tbl_size_ = 0;
}

size_t KhxReader::willneed_size() const {
    if (!mmap_.is_open()) return 0;
    return mmap_.size();
}

void KhxReader::apply_madvise(bool willneed) {
    if (!mmap_.is_open()) return;
    mmap_.advise(willneed ? MADV_WILLNEED : MADV_RANDOM);
}

uint64_t KhxReader::count_excluded() const {
    if (!is_open()) return 0;
    uint64_t count = 0;
    uint32_t bitset_bytes = (tbl_size_ + 7) / 8;
    for (uint32_t i = 0; i < bitset_bytes; i++) {
        count += __builtin_popcount(bitset_[i]);
    }
    return count;
}

} // namespace ikafssn
