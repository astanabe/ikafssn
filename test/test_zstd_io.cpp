// Round-trip test for src/util/zstd_oneshot.{hpp,cpp}.
//
// Exercises both the in-memory compress/decompress path and the
// `*_to_file` / `*_file` variants used by the async-job persistence
// (`<job_id>.deflines.zst`, `<job_id>.result.bin.zst`) and by the
// httpd ResultStore.

#include "util/zstd_oneshot.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <random>
#include <string>
#include <vector>

#include "test_util.hpp"

using ikafssn::zstd_compress;
using ikafssn::zstd_decompress;
using ikafssn::zstd_compress_to_file;
using ikafssn::zstd_decompress_file;

namespace fs = std::filesystem;

static std::string tmp_path() {
    auto dir = fs::temp_directory_path() / "ikafssn_zstd_io_test";
    std::error_code ec;
    fs::create_directories(dir, ec);
    static int counter = 0;
    return (dir / ("blob_" + std::to_string(++counter) + ".zst")).string();
}

static void test_round_trip_small() {
    const std::string original = "hello, ikafssn async world";
    std::vector<uint8_t> compressed;
    std::string err;
    CHECK(zstd_compress(original.data(), original.size(), compressed, 3, err));
    CHECK(!compressed.empty());

    std::vector<uint8_t> decompressed;
    CHECK(zstd_decompress(compressed.data(), compressed.size(),
                            decompressed, err));
    CHECK(decompressed.size() == original.size());
    CHECK(std::memcmp(decompressed.data(), original.data(),
                        original.size()) == 0);
}

static void test_round_trip_large_random() {
    std::mt19937 rng(42);
    std::uniform_int_distribution<int> dist(0, 255);
    std::vector<uint8_t> payload(256 * 1024);
    for (auto& b : payload) b = static_cast<uint8_t>(dist(rng));

    std::vector<uint8_t> compressed;
    std::string err;
    CHECK(zstd_compress(payload.data(), payload.size(),
                          compressed, 3, err));

    std::vector<uint8_t> decompressed;
    CHECK(zstd_decompress(compressed.data(), compressed.size(),
                            decompressed, err));
    CHECK(decompressed.size() == payload.size());
    CHECK(std::memcmp(decompressed.data(), payload.data(),
                        payload.size()) == 0);
}

static void test_round_trip_file_atomic() {
    const std::string deflines =
        "ACC0001.1 first defline\n"
        "ACC0002.1 second defline\n"
        "ACC0003.1 third defline\n";
    auto path = tmp_path();

    std::string err;
    CHECK(zstd_compress_to_file(path, deflines.data(), deflines.size(),
                                  3, err));
    CHECK(!fs::exists(path + ".tmp"));  // atomic rename cleaned up
    CHECK(fs::exists(path));

    std::vector<uint8_t> out;
    CHECK(zstd_decompress_file(path, out, err));
    CHECK(out.size() == deflines.size());
    CHECK(std::memcmp(out.data(), deflines.data(), deflines.size()) == 0);

    std::error_code ec;
    fs::remove(path, ec);
}

static void test_empty_payload() {
    std::vector<uint8_t> compressed;
    std::string err;
    CHECK(zstd_compress(nullptr, 0, compressed, 3, err));

    std::vector<uint8_t> decompressed;
    CHECK(zstd_decompress(compressed.data(), compressed.size(),
                            decompressed, err));
    CHECK(decompressed.empty());
}

int main() {
    test_round_trip_small();
    test_round_trip_large_random();
    test_round_trip_file_atomic();
    test_empty_payload();
    TEST_SUMMARY();
    return g_fail_count == 0 ? 0 : 1;
}
