#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace ikafssn {

// One-shot zstd compression/decompression helpers used by ikafssnclient
// for persisted defline lists and cached result blobs under
// `.ikafssnclient/<group_id>/`.  Both directions return false on failure
// and populate `error_msg`.

bool zstd_compress(const void* data, size_t size,
                   std::vector<uint8_t>& out,
                   int level,
                   std::string& error_msg);

bool zstd_decompress(const void* data, size_t size,
                     std::vector<uint8_t>& out,
                     std::string& error_msg);

bool zstd_compress_to_file(const std::string& path,
                           const void* data, size_t size,
                           int level,
                           std::string& error_msg);

bool zstd_decompress_file(const std::string& path,
                          std::vector<uint8_t>& out,
                          std::string& error_msg);

} // namespace ikafssn
