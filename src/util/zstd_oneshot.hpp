#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace ikafssn {

// One-shot zstd compression/decompression helpers shared by ikafssnclient
// (defline cache, cached result blobs) and ikafssnhttpd (per-job result
// files written by ResultStore, request-body decompression in HttpController).
// Both directions return false on failure and populate `error_msg`.

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
