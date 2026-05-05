#include "util/zstd_oneshot.hpp"

#include <cstdio>
#include <cstdlib>
#include <fstream>

#include <zstd.h>

namespace ikafssn {

bool zstd_compress(const void* data, size_t size,
                   std::vector<uint8_t>& out,
                   int level,
                   std::string& error_msg) {
    size_t bound = ZSTD_compressBound(size);
    out.resize(bound);
    size_t n = ZSTD_compress(out.data(), bound, data, size, level);
    if (ZSTD_isError(n)) {
        error_msg = "zstd compress failed: ";
        error_msg += ZSTD_getErrorName(n);
        return false;
    }
    out.resize(n);
    return true;
}

bool zstd_decompress(const void* data, size_t size,
                     std::vector<uint8_t>& out,
                     std::string& error_msg) {
    unsigned long long dec = ZSTD_getFrameContentSize(data, size);
    if (dec == ZSTD_CONTENTSIZE_ERROR) {
        error_msg = "zstd decompress: not a valid zstd frame";
        return false;
    }
    if (dec == ZSTD_CONTENTSIZE_UNKNOWN) {
        // Fallback: stream into a growing buffer.
        ZSTD_DStream* ds = ZSTD_createDStream();
        if (!ds) {
            error_msg = "zstd decompress: ZSTD_createDStream failed";
            return false;
        }
        ZSTD_initDStream(ds);
        ZSTD_inBuffer in{data, size, 0};
        out.clear();
        std::vector<uint8_t> chunk(64 * 1024);
        while (in.pos < in.size) {
            ZSTD_outBuffer ob{chunk.data(), chunk.size(), 0};
            size_t r = ZSTD_decompressStream(ds, &ob, &in);
            if (ZSTD_isError(r)) {
                error_msg = "zstd decompress stream failed: ";
                error_msg += ZSTD_getErrorName(r);
                ZSTD_freeDStream(ds);
                return false;
            }
            out.insert(out.end(), chunk.data(), chunk.data() + ob.pos);
            if (r == 0) break;
        }
        ZSTD_freeDStream(ds);
        return true;
    }
    out.resize(static_cast<size_t>(dec));
    size_t n = ZSTD_decompress(out.data(), out.size(), data, size);
    if (ZSTD_isError(n)) {
        error_msg = "zstd decompress failed: ";
        error_msg += ZSTD_getErrorName(n);
        return false;
    }
    out.resize(n);
    return true;
}

bool zstd_compress_to_file(const std::string& path,
                           const void* data, size_t size,
                           int level,
                           std::string& error_msg) {
    std::vector<uint8_t> blob;
    if (!zstd_compress(data, size, blob, level, error_msg)) return false;

    // Atomic write: tmp -> rename
    std::string tmp = path + ".tmp";
    {
        std::ofstream out(tmp, std::ios::binary | std::ios::trunc);
        if (!out.is_open()) {
            error_msg = "zstd_compress_to_file: cannot open " + tmp;
            return false;
        }
        out.write(reinterpret_cast<const char*>(blob.data()),
                  static_cast<std::streamsize>(blob.size()));
        if (!out.good()) {
            error_msg = "zstd_compress_to_file: write error to " + tmp;
            return false;
        }
    }
    if (std::rename(tmp.c_str(), path.c_str()) != 0) {
        error_msg = "zstd_compress_to_file: rename failed for " + path;
        std::remove(tmp.c_str());
        return false;
    }
    return true;
}

bool zstd_decompress_file(const std::string& path,
                          std::vector<uint8_t>& out,
                          std::string& error_msg) {
    std::ifstream in(path, std::ios::binary);
    if (!in.is_open()) {
        error_msg = "zstd_decompress_file: cannot open " + path;
        return false;
    }
    in.seekg(0, std::ios::end);
    auto sz = static_cast<size_t>(in.tellg());
    in.seekg(0, std::ios::beg);
    std::vector<uint8_t> blob(sz);
    if (sz > 0) {
        in.read(reinterpret_cast<char*>(blob.data()),
                static_cast<std::streamsize>(sz));
        if (!in.good() && !in.eof()) {
            error_msg = "zstd_decompress_file: read error from " + path;
            return false;
        }
    }
    return zstd_decompress(blob.data(), blob.size(), out, error_msg);
}

} // namespace ikafssn
