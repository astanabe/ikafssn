// Unit test for src/ikafssnhttpd/result_store.cpp.
//
// Covers: write/read round-trip with shard layout, atomic .tmp cleanup,
// sweep_tmp at startup, list_job_ids enumeration, unlink semantics.

#include "ikafssnhttpd/result_store.hpp"
#include "util/zstd_oneshot.hpp"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "test_util.hpp"

namespace fs = std::filesystem;

using ikafssn::ResultStore;

static std::string tmp_root() {
    auto dir = fs::temp_directory_path() / "ikafssn_result_store_test";
    std::error_code ec;
    fs::remove_all(dir, ec);
    fs::create_directories(dir, ec);
    return dir.string();
}

static void test_write_read_round_trip() {
    auto root = tmp_root();
    ResultStore rs(root, 3);
    std::string err;
    CHECK(rs.init(err));

    std::string job_id = "abcd1234-5678-4abc-89ab-cdef01234567";
    std::string payload(8192, '\0');
    for (size_t i = 0; i < payload.size(); ++i) {
        payload[i] = static_cast<char>(i & 0xFF);
    }

    CHECK(rs.write(job_id, payload.data(), payload.size(), err));
    CHECK(rs.exists(job_id));

    std::string path = rs.path_for(job_id);
    // 2-char prefix shard.
    std::string expected = root + "/ab/" + job_id + ".bin.zst";
    CHECK(path == expected);
    CHECK(!fs::exists(path + ".tmp"));

    // Decompress directly to verify the file is a single zstd frame.
    std::vector<uint8_t> decoded;
    CHECK(ikafssn::zstd_decompress_file(path, decoded, err));
    CHECK(decoded.size() == payload.size());
    CHECK(std::memcmp(decoded.data(), payload.data(), payload.size()) == 0);
}

static void test_overwrite() {
    auto root = tmp_root();
    ResultStore rs(root, 3);
    std::string err;
    CHECK(rs.init(err));

    std::string job_id = "ff000000-0000-4000-8000-000000000001";
    CHECK(rs.write(job_id, "first", 5, err));
    CHECK(rs.write(job_id, "second-overwrite", 16, err));

    std::vector<uint8_t> decoded;
    CHECK(ikafssn::zstd_decompress_file(rs.path_for(job_id), decoded, err));
    CHECK(decoded.size() == 16);
    CHECK(std::memcmp(decoded.data(), "second-overwrite", 16) == 0);
}

static void test_unlink_missing_is_ok() {
    auto root = tmp_root();
    ResultStore rs(root, 3);
    std::string err;
    CHECK(rs.init(err));

    std::string job_id = "ee000000-0000-4000-8000-000000000002";
    // Never written; unlink should still succeed.
    CHECK(rs.unlink(job_id, err));
}

static void test_sweep_tmp() {
    auto root = tmp_root();
    {
        // Pre-create a stranded .tmp file that init() should sweep.
        fs::create_directories(root + "/aa");
        std::ofstream f(root + "/aa/garbage.bin.zst.tmp", std::ios::binary);
        f << "leftover from a crash";
    }
    ResultStore rs(root, 3);
    std::string err;
    CHECK(rs.init(err));
    CHECK(!fs::exists(root + "/aa/garbage.bin.zst.tmp"));
}

static void test_list_job_ids() {
    auto root = tmp_root();
    ResultStore rs(root, 3);
    std::string err;
    CHECK(rs.init(err));

    std::vector<std::string> ids = {
        "11111111-aaaa-4000-8000-000000000001",
        "22222222-bbbb-4000-8000-000000000002",
        "33333333-cccc-4000-8000-000000000003",
    };
    for (const auto& jid : ids) {
        CHECK(rs.write(jid, jid.data(), jid.size(), err));
    }

    auto found = rs.list_job_ids(err);
    CHECK(err.empty());
    std::sort(found.begin(), found.end());
    auto expected = ids;
    std::sort(expected.begin(), expected.end());
    CHECK(found.size() == expected.size());
    for (size_t i = 0; i < expected.size(); ++i) {
        CHECK(found[i] == expected[i]);
    }
}

int main() {
    test_write_read_round_trip();
    test_overwrite();
    test_unlink_missing_is_ok();
    test_sweep_tmp();
    test_list_job_ids();
    TEST_SUMMARY();
    return g_fail_count == 0 ? 0 : 1;
}
