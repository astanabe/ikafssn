#include "index/kix_reader.hpp"
#include "core/config.hpp"
#include "core/varint.hpp"

#include <sys/mman.h>
#include <cstring>
#include <cstdio>
#include <vector>

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

    offset32_ = (header_->flags & KIX_FLAG_OFFSET32) != 0;

    const uint8_t* ptr = mmap_.data() + sizeof(KixHeader);

    // offsets has table_size_ + 1 entries (sentinel at end)
    if (offset32_) {
        offsets32_ = reinterpret_cast<const uint32_t*>(ptr);
        ptr += sizeof(uint32_t) * (table_size_ + 1);
    } else {
        offsets64_ = reinterpret_cast<const uint64_t*>(ptr);
        ptr += sizeof(uint64_t) * (table_size_ + 1);
    }

    posting_data_ = ptr;
    posting_data_size_ = mmap_.size() - (ptr - mmap_.data());

    return true;
}

void KixReader::close() {
    mmap_.close();
    header_ = nullptr;
    offsets64_ = nullptr;
    offsets32_ = nullptr;
    offset32_ = false;
    posting_data_ = nullptr;
    posting_data_size_ = 0;
    table_size_ = 0;
}

size_t KixReader::willneed_size() const {
    if (!mmap_.is_open()) return 0;
    size_t offset_bytes = offset32_
        ? sizeof(uint32_t) * (table_size_ + 1)
        : sizeof(uint64_t) * (table_size_ + 1);
    return sizeof(KixHeader) + offset_bytes;
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

uint32_t KixReader::count_postings(uint32_t kmer) const {
    uint64_t byte_len = posting_byte_length(kmer);
    if (byte_len == 0) return 0;
    const uint8_t* ptr = posting_data_ + posting_offset(kmer);
    const uint8_t* end = ptr + byte_len;
    uint32_t count = 0;
    while (ptr < end) {
        uint32_t dummy;
        ptr += varint_decode(ptr, dummy);
        count++;
    }
    return count;
}

std::vector<uint32_t> KixReader::bulk_count_postings() const {
    std::vector<uint32_t> counts(table_size_, 0);
    for (uint32_t i = 0; i < table_size_; i++) {
        counts[i] = count_postings(i);
    }
    return counts;
}

} // namespace ikafssn
