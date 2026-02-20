#include "index/ksx_reader.hpp"
#include "index/ksx_format.hpp"
#include "core/config.hpp"

#include <sys/mman.h>
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

    if (std::memcmp(hdr->magic, KSX_MAGIC, 4) != 0) {
        std::fprintf(stderr, "KsxReader: invalid magic\n");
        close();
        return false;
    }

    if (hdr->format_version != KSX_FORMAT_VERSION) {
        std::fprintf(stderr, "KsxReader: unsupported format version %u\n", hdr->format_version);
        close();
        return false;
    }

    num_sequences_ = hdr->num_sequences;

    const uint8_t* ptr = mmap_.data() + sizeof(KsxHeader);
    seq_lengths_ = reinterpret_cast<const uint32_t*>(ptr);
    ptr += sizeof(uint32_t) * num_sequences_;

    acc_offsets_ = reinterpret_cast<const uint32_t*>(ptr);
    ptr += sizeof(uint32_t) * (num_sequences_ + 1);

    acc_strings_ = reinterpret_cast<const char*>(ptr);

    // Accession lookups are random (by OID from search results)
    mmap_.advise(MADV_RANDOM);

    return true;
}

void KsxReader::close() {
    mmap_.close();
    num_sequences_ = 0;
    seq_lengths_ = nullptr;
    acc_offsets_ = nullptr;
    acc_strings_ = nullptr;
}

uint32_t KsxReader::seq_length(uint32_t oid) const {
    return seq_lengths_[oid];
}

std::string_view KsxReader::accession(uint32_t oid) const {
    uint32_t start = acc_offsets_[oid];
    uint32_t len = acc_offsets_[oid + 1] - start;
    return std::string_view(acc_strings_ + start, len);
}

} // namespace ikafssn
