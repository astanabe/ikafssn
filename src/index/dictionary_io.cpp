#include "index/dictionary_io.hpp"
#include "index/ef_codec.hpp"

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <vector>

namespace ikafssn {

namespace {

bool write_blob(std::FILE* fp, const std::uint64_t* offsets, std::size_t D,
                std::uint64_t posting_file_size) {
    std::vector<std::uint8_t> blob;
    ef::encode_dictionary_ef(offsets, D, posting_file_size, blob);
    if (blob.empty()) return false;
    return std::fwrite(blob.data(), 1, blob.size(), fp) == blob.size();
}

} // anonymous namespace

bool write_kix_dictionary_ef(std::FILE* fp,
                             const std::uint64_t* offsets,
                             std::uint32_t table_size,
                             std::uint64_t posting_file_size) {
    return write_blob(fp, offsets,
                      static_cast<std::size_t>(table_size) + 1,
                      posting_file_size);
}

bool write_kpx_dictionary_ef(std::FILE* fp,
                             const std::uint64_t* pos_offsets,
                             std::uint32_t table_size,
                             std::uint64_t posting_file_size) {
    return write_blob(fp, pos_offsets,
                      static_cast<std::size_t>(table_size),
                      posting_file_size);
}

} // namespace ikafssn
