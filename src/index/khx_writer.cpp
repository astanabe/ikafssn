#include "index/khx_writer.hpp"
#include "index/khx_format.hpp"
#include "core/config.hpp"
#include "util/logger.hpp"

#include <cstdio>
#include <cstring>

namespace ikafssn {

bool write_khx(const std::string& path, int k,
               const std::vector<uint32_t>& counts,
               uint64_t freq_threshold,
               const Logger& logger) {

    uint64_t tbl_size = table_size(k);
    if (counts.size() != tbl_size) {
        logger.error("write_khx: counts size mismatch (got %zu, expected %lu)",
                     counts.size(), static_cast<unsigned long>(tbl_size));
        return false;
    }

    FILE* fp = std::fopen(path.c_str(), "wb");
    if (!fp) {
        logger.error("write_khx: cannot open %s for writing", path.c_str());
        return false;
    }

    // Write header
    KhxHeader hdr{};
    std::memcpy(hdr.magic, KHX_MAGIC, 4);
    hdr.format_version = KHX_FORMAT_VERSION;
    hdr.k = static_cast<uint8_t>(k);
    std::fwrite(&hdr, sizeof(hdr), 1, fp);

    // Build bitset: ceil(tbl_size / 8) bytes
    uint64_t bitset_bytes = (tbl_size + 7) / 8;
    std::vector<uint8_t> bitset(bitset_bytes, 0);

    uint64_t excluded_count = 0;
    for (uint64_t i = 0; i < tbl_size; i++) {
        if (counts[i] > freq_threshold) {
            bitset[i / 8] |= (1u << (i % 8));
            excluded_count++;
        }
    }

    std::fwrite(bitset.data(), 1, bitset_bytes, fp);
    std::fclose(fp);

    logger.info("Wrote %s (%lu excluded k-mers)", path.c_str(),
                static_cast<unsigned long>(excluded_count));
    return true;
}

} // namespace ikafssn
