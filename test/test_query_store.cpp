// Unit test for src/ikafssnhttpd/query_store.cpp.
//
// Covers:
//   - write_compressed_passthrough preserves the on-disk bytes verbatim
//     (i.e. an externally-produced zstd frame with a different
//     compression level survives round-trip without re-compression);
//   - write_plain + read round-trips a plaintext JSON body through the
//     configured zstd level;
//   - 2-character shard layout + path_for;
//   - sweep_tmp removes stranded *.json.zst.tmp files at init();
//   - list_job_ids enumerates persisted job IDs and skips .tmp files;
//   - unlink on a missing file is a no-op success.

#include "ikafssnhttpd/query_store.hpp"
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

using ikafssn::QueryStore;

static std::string tmp_root() {
    auto dir = fs::temp_directory_path() / "ikafssn_query_store_test";
    std::error_code ec;
    fs::remove_all(dir, ec);
    fs::create_directories(dir, ec);
    return dir.string();
}

static void test_passthrough_round_trip() {
    auto root = tmp_root();
    QueryStore qs(root, /*level=*/3);
    std::string err;
    CHECK(qs.init(err));

    std::string job_id = "abcd1234-5678-4abc-89ab-cdef01234567";
    std::string body =
        "{\"job_id\":\"" + job_id + "\",\"queries\":[{\"qseqid\":\"q1\","
        "\"sequence\":\"ACGT\"}]}";

    // Externally compress with level 19 — store passthrough should
    // leave the bytes alone, not recompress at the QueryStore level.
    std::vector<uint8_t> precompressed;
    CHECK(ikafssn::zstd_compress(body.data(), body.size(),
                                  precompressed, 19, err));

    CHECK(qs.write_compressed_passthrough(
        job_id, precompressed.data(), precompressed.size(), err));
    CHECK(qs.exists(job_id));

    std::string path = qs.path_for(job_id);
    std::string expected = root + "/ab/" + job_id + ".json.zst";
    CHECK(path == expected);
    CHECK(!fs::exists(path + ".tmp"));

    // The bytes on disk must be identical to what we wrote.
    std::ifstream f(path, std::ios::binary);
    std::vector<uint8_t> on_disk((std::istreambuf_iterator<char>(f)),
                                  std::istreambuf_iterator<char>());
    CHECK(on_disk.size() == precompressed.size());
    CHECK(std::memcmp(on_disk.data(), precompressed.data(),
                      precompressed.size()) == 0);

    // QueryStore::read should decompress it back to the original JSON.
    std::string round_trip;
    CHECK(qs.read(job_id, round_trip, err));
    CHECK(round_trip == body);
}

static void test_plain_round_trip() {
    auto root = tmp_root();
    QueryStore qs(root, /*level=*/3);
    std::string err;
    CHECK(qs.init(err));

    std::string job_id = "11111111-aaaa-4000-8000-000000000001";
    std::string body =
        "{\"job_id\":\"" + job_id + "\",\"queries\":[{\"qseqid\":\"q\","
        "\"sequence\":\"ACGTACGT\"}]}";

    CHECK(qs.write_plain(job_id, body.data(), body.size(), err));
    CHECK(qs.exists(job_id));

    std::string round_trip;
    CHECK(qs.read(job_id, round_trip, err));
    CHECK(round_trip == body);

    // path_for / shard_for sanity.
    CHECK(qs.path_for(job_id) ==
          root + "/11/" + job_id + ".json.zst");
}

static void test_unlink_missing_is_ok() {
    auto root = tmp_root();
    QueryStore qs(root, 3);
    std::string err;
    CHECK(qs.init(err));

    std::string job_id = "ee000000-0000-4000-8000-000000000002";
    CHECK(qs.unlink(job_id, err));
}

static void test_sweep_tmp() {
    auto root = tmp_root();
    {
        // Pre-create a stranded .tmp file that init() should sweep.
        fs::create_directories(root + "/aa");
        std::ofstream f(root + "/aa/garbage.json.zst.tmp", std::ios::binary);
        f << "leftover from a crash";
    }
    QueryStore qs(root, 3);
    std::string err;
    CHECK(qs.init(err));
    CHECK(!fs::exists(root + "/aa/garbage.json.zst.tmp"));
}

static void test_list_job_ids() {
    auto root = tmp_root();
    QueryStore qs(root, 3);
    std::string err;
    CHECK(qs.init(err));

    std::vector<std::string> ids = {
        "11111111-aaaa-4000-8000-000000000001",
        "22222222-bbbb-4000-8000-000000000002",
        "33333333-cccc-4000-8000-000000000003",
    };
    std::string body = "{\"job_id\":\"x\",\"queries\":[{\"qseqid\":\"q\","
                       "\"sequence\":\"A\"}]}";
    for (const auto& jid : ids) {
        CHECK(qs.write_plain(jid, body.data(), body.size(), err));
    }
    // Drop a stranded *.json.zst.tmp; list_job_ids must skip it.
    fs::create_directories(root + "/de");
    {
        std::ofstream f(root + "/de/should-be-skipped.json.zst.tmp",
                        std::ios::binary);
        f << "garbage";
    }

    auto found = qs.list_job_ids(err);
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
    test_passthrough_round_trip();
    test_plain_round_trip();
    test_unlink_missing_is_ok();
    test_sweep_tmp();
    test_list_job_ids();
    TEST_SUMMARY();
    return g_fail_count == 0 ? 0 : 1;
}
