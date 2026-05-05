#include "index/kix_reader.hpp"
#include "core/config.hpp"
#include "index/pfd_codec.hpp"

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
        std::fprintf(stderr,
            "KixReader: index format version mismatch (got %u, expected %u). "
            "Please rebuild with the current ikafssnindex.\n",
            header_->format_version, KIX_FORMAT_VERSION);
        close();
        return false;
    }

    table_size_ = ikafssn::table_size(header_->k);

    // Phase 7a: Elias-Fano dictionary follows the header.  The legacy
    // KIX_FLAG_OFFSET32 bit is ignored at read time.
    const uint8_t* dict_ptr = mmap_.data() + sizeof(KixHeader);
    const size_t   remaining = mmap_.size() - sizeof(KixHeader);
    if (!dict_.open(dict_ptr, remaining)) {
        std::fprintf(stderr, "KixReader: invalid Elias-Fano dictionary\n");
        close();
        return false;
    }
    if (dict_.size() != static_cast<size_t>(table_size_) + 1) {
        std::fprintf(stderr,
            "KixReader: dictionary entry count %zu does not match expected %zu\n",
            dict_.size(), static_cast<size_t>(table_size_) + 1);
        close();
        return false;
    }

    const uint8_t* ptr = dict_ptr + dict_.blob_bytes();
    posting_file_ = ptr;
    posting_file_size_ = mmap_.size() - (ptr - mmap_.data());

    return true;
}

void KixReader::close() {
    mmap_.close();
    header_ = nullptr;
    dict_ = ef::EFDictionary{};
    posting_file_ = nullptr;
    posting_file_size_ = 0;
    table_size_ = 0;
}

size_t KixReader::willneed_size() const {
    if (!mmap_.is_open()) return 0;
    return sizeof(KixHeader) + dict_.willneed_size();
}

void KixReader::apply_madvise(bool willneed) {
    mmap_.advise_dict_posting(willneed_size(), willneed);
}

size_t KixReader::willneed_size_full() const {
    if (!mmap_.is_open()) return 0;
    return mmap_.size();
}

void KixReader::apply_madvise_full(bool willneed) {
    if (!mmap_.is_open()) return;
    const size_t dict_size = willneed_size();
    if (willneed) {
        mmap_.advise(0, dict_size, MADV_WILLNEED);
#ifdef MADV_HUGEPAGE
        mmap_.advise(0, dict_size, MADV_HUGEPAGE);
#endif
        if (mmap_.size() > dict_size)
            mmap_.advise(dict_size, mmap_.size() - dict_size, MADV_WILLNEED);
    } else {
        if (mmap_.size() > dict_size)
            mmap_.advise(dict_size, mmap_.size() - dict_size, MADV_DONTNEED);
    }
}

void KixReader::apply_madvise_posting_random() {
    if (!mmap_.is_open()) return;
    const size_t dict_size = willneed_size();
    mmap_.advise(0, dict_size, MADV_WILLNEED);
#ifdef MADV_HUGEPAGE
    mmap_.advise(0, dict_size, MADV_HUGEPAGE);
#endif
    if (mmap_.size() > dict_size)
        mmap_.advise(dict_size, mmap_.size() - dict_size, MADV_RANDOM);
}

std::vector<uint32_t> KixReader::bulk_count_postings() const {
    std::vector<uint32_t> counts(table_size_, 0);
    // Walk the dictionary once with sequential access_pair calls so the
    // EF select1 sample lookup is amortised across consecutive k-mers.
    // (7b SIMD will turn this into a streaming decode; for the 7a scalar
    // PoC the per-call cost is acceptable since this path runs only at
    // ikafssnindex finalize / ikafssninfo time.)
    uint64_t s = dict_.access(0);
    for (uint32_t i = 0; i < table_size_; i++) {
        uint64_t e = dict_.access(i + 1);
        uint64_t byte_len = e - s;
        if (byte_len > 0) {
            counts[i] = pfd::posting_count(posting_file_ + s, byte_len);
        }
        s = e;
    }
    return counts;
}

} // namespace ikafssn
