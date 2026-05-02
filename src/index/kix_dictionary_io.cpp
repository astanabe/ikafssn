#include "index/kix_dictionary_io.hpp"
#include "index/ef_codec.hpp"

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <vector>

namespace ikafssn {

bool write_kix_dictionary_ef(std::FILE* fp,
                             const std::uint64_t* offsets,
                             std::uint32_t table_size,
                             std::uint64_t posting_file_size) {
    std::vector<std::uint8_t> blob;
    const std::size_t D = static_cast<std::size_t>(table_size) + 1;
    ef::encode_dictionary_ef(offsets, D, posting_file_size, blob);
    if (blob.empty()) {
        return false;
    }
    if (std::fwrite(blob.data(), 1, blob.size(), fp) != blob.size()) {
        return false;
    }
    return true;
}

} // namespace ikafssn
