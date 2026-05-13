#include "index/volume_validator.hpp"

#include "index/index_filter.hpp"
#include "index/khx_format.hpp"
#include "index/kix_format.hpp"
#include "index/kix_reader.hpp"
#include "index/ksx_format.hpp"
#include "core/config.hpp"
#include "util/logger.hpp"

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <sys/stat.h>

namespace ikafssn {

namespace {

int64_t file_size_or_neg(const std::string& path) {
    struct stat st;
    if (::stat(path.c_str(), &st) != 0) return -1;
    return static_cast<int64_t>(st.st_size);
}

bool verify_ksx_layout(const std::string& ksx_path,
                       const Logger& logger) {
    FILE* fp = std::fopen(ksx_path.c_str(), "rb");
    if (!fp) return false;
    KsxHeader hdr{};
    if (std::fread(&hdr, sizeof(hdr), 1, fp) != 1) {
        std::fclose(fp);
        return false;
    }
    if (std::memcmp(hdr.magic, KSX_MAGIC, sizeof(KSX_MAGIC)) != 0
        || hdr.format_version != KSX_FORMAT_VERSION) {
        std::fclose(fp);
        logger.warn("validate: %s magic / format_version mismatch", ksx_path.c_str());
        return false;
    }
    const uint32_t num_sequences = hdr.num_sequences;
    const uint32_t num_parents   = hdr.num_parents;

    const long acc_offsets_pos = static_cast<long>(sizeof(KsxHeader))
                               + static_cast<long>(sizeof(uint32_t)) * num_parents
                               + static_cast<long>(sizeof(uint32_t)) * num_parents
                               + static_cast<long>(sizeof(uint32_t)) * num_parents;
    if (std::fseek(fp, acc_offsets_pos, SEEK_SET) != 0) {
        std::fclose(fp);
        return false;
    }
    uint32_t string_table_bytes = 0;
    if (std::fread(&string_table_bytes, sizeof(uint32_t), 1, fp) != 1) {
        std::fclose(fp);
        return false;
    }
    std::fclose(fp);

    int64_t expected = static_cast<int64_t>(sizeof(KsxHeader))
                     + static_cast<int64_t>(sizeof(uint32_t)) * num_parents
                     + static_cast<int64_t>(sizeof(uint32_t)) * num_parents
                     + static_cast<int64_t>(sizeof(uint32_t)) * (num_parents + 1)
                     + static_cast<int64_t>(string_table_bytes)
                     + static_cast<int64_t>(sizeof(uint32_t)) * num_sequences
                     + static_cast<int64_t>(sizeof(uint32_t)) * num_sequences
                     + static_cast<int64_t>(sizeof(uint32_t)) * num_sequences;
    int64_t actual = file_size_or_neg(ksx_path);
    if (actual != expected) {
        logger.warn("validate: %s file size %ld != expected %ld",
                    ksx_path.c_str(),
                    static_cast<long>(actual),
                    static_cast<long>(expected));
        return false;
    }
    return true;
}

// Confirms the posting-file region matches the EF dictionary's
// sentinel offset (catches both truncation and oversized files).
bool verify_kix_layout(const std::string& kix_path,
                       const Logger& logger) {
    KixReader kix;
    if (!kix.open(kix_path)) {
        logger.warn("validate: cannot open %s", kix_path.c_str());
        return false;
    }
    uint64_t expected_post = kix.posting_list_offset(kix.table_size());
    uint64_t actual_post = kix.posting_file_size();
    if (expected_post != actual_post) {
        logger.warn("validate: %s posting file size %lu != EF sentinel %lu",
                    kix_path.c_str(),
                    static_cast<unsigned long>(actual_post),
                    static_cast<unsigned long>(expected_post));
        return false;
    }
    return true;
}

bool validate_volume_strict_impl(const std::string& prefix,
                                 const char* ksx_suffix,
                                 const char* kix_suffix,
                                 const char* kpx_suffix,
                                 bool skip_kpx,
                                 const Logger& logger) {
    const std::string ksx_path = prefix + ksx_suffix;
    const std::string kix_path = prefix + kix_suffix;
    const std::string kpx_path = skip_kpx
        ? std::string{} : (prefix + kpx_suffix);

    if (!std::filesystem::exists(ksx_path)) return false;
    if (!std::filesystem::exists(kix_path)) return false;
    if (!skip_kpx && !std::filesystem::exists(kpx_path)) return false;

    if (!validate_ksx_file_strict(ksx_path, logger)) return false;
    if (!validate_kix_kpx_strict(kix_path, kpx_path, logger)) return false;
    return true;
}

} // anonymous namespace

bool validate_ksx_file_strict(const std::string& ksx_path,
                              const Logger& logger) {
    return verify_ksx_layout(ksx_path, logger);
}

bool validate_kix_kpx_strict(const std::string& kix_path,
                             const std::string& kpx_path,
                             const Logger& logger) {
    if (!verify_kix_layout(kix_path, logger)) return false;
    if (!validate_volume(kix_path, kpx_path, nullptr, logger)) {
        return false;
    }
    return true;
}

bool validate_khx_file_strict(const std::string& khx_path,
                              int k, uint8_t t, uint8_t template_type,
                              const Logger& logger) {
    FILE* fp = std::fopen(khx_path.c_str(), "rb");
    if (!fp) return false;

    KhxHeader hdr{};
    bool ok = (std::fread(&hdr, sizeof(hdr), 1, fp) == 1);
    std::fclose(fp);
    if (!ok) return false;

    if (std::memcmp(hdr.magic, KHX_MAGIC, 4) != 0) {
        logger.warn("validate: %s has invalid magic", khx_path.c_str());
        return false;
    }
    if (hdr.format_version != KHX_FORMAT_VERSION) {
        logger.warn("validate: %s format_version=%u != expected %u",
                    khx_path.c_str(), hdr.format_version, KHX_FORMAT_VERSION);
        return false;
    }
    if (hdr.k != static_cast<uint8_t>(k) ||
        hdr.t != t ||
        hdr.template_type != template_type) {
        logger.warn("validate: %s (k=%u, t=%u, template_type=%u) does not match "
                    "requested (k=%d, t=%u, template_type=%u)",
                    khx_path.c_str(), hdr.k, hdr.t, hdr.template_type,
                    k, t, template_type);
        return false;
    }

    int64_t expected = static_cast<int64_t>(sizeof(KhxHeader))
                     + static_cast<int64_t>((table_size(k) + 7) / 8);
    int64_t actual = file_size_or_neg(khx_path);
    if (actual != expected) {
        logger.warn("validate: %s file size %ld != expected %ld",
                    khx_path.c_str(),
                    static_cast<long>(actual),
                    static_cast<long>(expected));
        return false;
    }
    return true;
}

bool validate_volume_final_strict(const std::string& prefix,
                                  bool skip_kpx,
                                  const Logger& logger) {
    return validate_volume_strict_impl(prefix, ".ksx", ".kix", ".kpx",
                                       skip_kpx, logger);
}

bool validate_volume_tmp_strict(const std::string& prefix,
                                bool skip_kpx,
                                const Logger& logger) {
    return validate_volume_strict_impl(prefix, ".ksx.tmp", ".kix.tmp",
                                       ".kpx.tmp", skip_kpx, logger);
}

} // namespace ikafssn
