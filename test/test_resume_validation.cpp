#include "test_util.hpp"
#include "ssu_test_fixture.hpp"
#include "index/volume_validator.hpp"
#include "io/volume_discovery.hpp"
#include "util/logger.hpp"

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>
#include <sys/wait.h>

using namespace ikafssn;
using namespace ssu_fixture;

namespace {

constexpr int kTestK = 7;
constexpr uint32_t kMinSeqLength = 64;

std::string g_testdb_path;
std::string g_index_dir;
std::string g_build_dir;
std::string g_vol_basename;

// Final-file stem for the single SSU volume; with the test flags
// (-max_freq_build defaulted, -max_degen_expand 0) no maxfreq /
// maxexpand suffix is appended.
std::string final_prefix() {
    return index_file_stem(g_index_dir, g_vol_basename, kTestK,
                           /*t=*/0, /*template_type=*/0,
                           kMinSeqLength,
                           /*min_length_split=*/50000,
                           /*overlap_length=*/500,
                           /*max_freq_build=*/1,
                           /*max_degen_expand=*/0);
}

void reset_index_dir() {
    std::error_code ec;
    std::filesystem::remove_all(g_index_dir, ec);
    std::filesystem::create_directories(g_index_dir);
}

int run_ikafssnindex(const std::string& extra_args, std::string* captured = nullptr) {
    std::string cmd = g_build_dir + "/ikafssnindex"
        + " -db " + g_testdb_path
        + " -k " + std::to_string(kTestK)
        + " -o " + g_index_dir
        + " -mode 1 -max_degen_expand 0 "
        + extra_args
        + " 2>&1";
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) return -1;
    char buf[4096];
    std::string out;
    while (std::fgets(buf, sizeof(buf), pipe)) {
        out += buf;
    }
    int ret = pclose(pipe);
    if (captured) *captured = out;
    if (WIFEXITED(ret)) return WEXITSTATUS(ret);
    return -1;
}

// XOR one byte at the given offset in-place on disk.
void flip_byte(const std::string& path, std::streamoff off, char xor_mask) {
    std::fstream fs(path, std::ios::in | std::ios::out | std::ios::binary);
    CHECK(fs.is_open());
    fs.seekg(off);
    char b = 0;
    fs.read(&b, 1);
    b ^= xor_mask;
    fs.seekp(off);
    fs.write(&b, 1);
    fs.close();
}

// Number of posting file scratch files (.kix.pf.tmp / .kpx.pf.tmp)
// currently in the index directory.
int count_pf_tmp() {
    int n = 0;
    for (auto& e : std::filesystem::directory_iterator(g_index_dir)) {
        const std::string name = e.path().filename().string();
        if (name.size() >= 7 && name.compare(name.size() - 7, 7, ".pf.tmp") == 0) {
            n++;
        }
    }
    return n;
}

void write_garbage(const std::string& path) {
    std::ofstream fs(path, std::ios::binary | std::ios::trunc);
    CHECK(fs.is_open());
    const std::string junk(4096, '\xA5');
    fs.write(junk.data(), static_cast<std::streamsize>(junk.size()));
    fs.close();
}

void setup() {
    check_ssu_available();

    g_testdb_path = ssu_db_prefix();
    g_index_dir = test_tmpdir("/tmp/ikafssn_resume_test");
    g_build_dir = std::string(SOURCE_DIR) + "/build/src";
    g_vol_basename = std::filesystem::path(g_testdb_path).filename().string();

    if (!std::filesystem::exists(g_build_dir + "/ikafssnindex")) {
        skip("ikafssnindex binary not found (expected under build/src/)");
    }

    reset_index_dir();
}

// Re-running on a pristine index reuses every volume.
void test_resume_full_complete_skip() {
    std::fprintf(stderr, "-- test_resume_full_complete_skip\n");
    reset_index_dir();

    CHECK_EQ(run_ikafssnindex(""), 0);

    std::vector<std::filesystem::file_time_type> before;
    for (auto& entry : std::filesystem::directory_iterator(g_index_dir)) {
        before.push_back(entry.last_write_time());
    }

    std::string out;
    CHECK_EQ(run_ikafssnindex("", &out), 0);
    CHECK(out.find("=== Collecting metadata") == std::string::npos);
    CHECK(out.find("=== Writing postings") == std::string::npos);

    // The .kvx manifest is re-written every run; exclude it from the
    // mtime comparison.
    int unchanged_count = 0;
    for (auto& entry : std::filesystem::directory_iterator(g_index_dir)) {
        std::string ext = entry.path().extension().string();
        if (ext == ".kvx") continue;
        bool found = false;
        for (auto t : before) {
            if (t == entry.last_write_time()) { found = true; break; }
        }
        CHECK(found);
        if (found) unchanged_count++;
    }
    CHECK(unchanged_count > 0);
}

// Missing .kix triggers a per-volume rebuild.
void test_resume_one_volume_missing() {
    std::fprintf(stderr, "-- test_resume_one_volume_missing\n");
    reset_index_dir();

    CHECK_EQ(run_ikafssnindex(""), 0);
    std::filesystem::remove(final_prefix() + ".kix");

    std::string out;
    CHECK_EQ(run_ikafssnindex("", &out), 0);
    CHECK(out.find("=== Collecting metadata") != std::string::npos);
    CHECK(out.find("=== Writing postings") != std::string::npos);
    CHECK(std::filesystem::exists(final_prefix() + ".kix"));
}

// Corrupted .kix is rejected by strict validation and rebuilt.
void test_resume_one_volume_corrupted() {
    std::fprintf(stderr, "-- test_resume_one_volume_corrupted\n");
    reset_index_dir();

    CHECK_EQ(run_ikafssnindex(""), 0);

    // Flip one byte at the start of the EF dictionary (just past the
    // 96 B .kix header).
    std::string kix = final_prefix() + ".kix";
    flip_byte(kix, /*off=*/96, /*xor_mask=*/0x10);

    std::string out;
    CHECK_EQ(run_ikafssnindex("", &out), 0);
    CHECK(out.find("failed strict validation") != std::string::npos);
    CHECK(out.find("=== Collecting metadata") != std::string::npos);
    CHECK(out.find("=== Writing postings") != std::string::npos);

    Logger silent(Logger::kError);
    CHECK(validate_volume_final_strict(final_prefix(),
                                       /*skip_kpx=*/true, silent));
}

// A volume left with only .ksx.tmp (metadata done, postings pending)
// resumes with the postings pass alone.
void test_resume_metadata_tmp_only() {
    std::fprintf(stderr, "-- test_resume_metadata_tmp_only\n");
    reset_index_dir();

    CHECK_EQ(run_ikafssnindex(""), 0);

    // Simulate a crash after build_metadata but before build_postings.
    std::string fp = final_prefix();
    std::filesystem::rename(fp + ".ksx", fp + ".ksx.tmp");
    std::filesystem::remove(fp + ".kix");

    std::string out;
    CHECK_EQ(run_ikafssnindex("", &out), 0);
    CHECK(out.find("reuse-ksx=1") != std::string::npos);
    CHECK(out.find("=== Writing postings") != std::string::npos);
}

// -force_rebuild 1 deletes and rewrites everything.
void test_force_rebuild_overrides() {
    std::fprintf(stderr, "-- test_force_rebuild_overrides\n");
    reset_index_dir();

    CHECK_EQ(run_ikafssnindex(""), 0);

    auto kix_path = final_prefix() + ".kix";
    auto before = std::filesystem::last_write_time(kix_path);

    std::string out;
    CHECK_EQ(run_ikafssnindex("-force_rebuild 1", &out), 0);
    auto after = std::filesystem::last_write_time(kix_path);
    CHECK(after != before);
}

// Fractional -max_freq_build resume: the file name encodes the
// resolved absolute threshold, and the peek pass must recover it from
// the existing .ksx so the fast-path skips metadata / postings / filter.
void test_resume_fractional_freq_filter() {
    std::fprintf(stderr, "-- test_resume_fractional_freq_filter\n");
    reset_index_dir();

    CHECK_EQ(run_ikafssnindex("-max_freq_build 0.01"), 0);

    std::string out;
    CHECK_EQ(run_ikafssnindex("-max_freq_build 0.01", &out), 0);
    CHECK(out.find("Resolved -max_freq_build=0.01") != std::string::npos);
    CHECK(out.find("=== Collecting metadata") == std::string::npos);
    CHECK(out.find("=== Writing postings") == std::string::npos);
    CHECK(out.find("Cross-volume filter:") == std::string::npos);
}

// .khx alone deleted with freq_filter_active: every volume must be
// demoted to kNone so the filter pass can regenerate .khx from fresh
// .kix.tmp inputs.
void test_resume_khx_only_invalid() {
    std::fprintf(stderr, "-- test_resume_khx_only_invalid\n");
    reset_index_dir();

    CHECK_EQ(run_ikafssnindex("-max_freq_build 100"), 0);

    std::string khx_path;
    for (auto& e : std::filesystem::directory_iterator(g_index_dir)) {
        if (e.path().extension() == ".khx") {
            khx_path = e.path().string();
            break;
        }
    }
    CHECK(!khx_path.empty());
    std::filesystem::remove(khx_path);

    std::string out;
    CHECK_EQ(run_ikafssnindex("-max_freq_build 100", &out), 0);
    CHECK(out.find("failed strict validation") != std::string::npos);
    CHECK(out.find("Demoting volume") != std::string::npos);
    CHECK(out.find("Cross-volume filter:") != std::string::npos);
    CHECK(std::filesystem::exists(khx_path));
}

// Direct API test: a corrupted .kix is rejected by
// validate_volume_final_strict (the post-build validation entry point).
void test_post_build_validation_failure_api() {
    std::fprintf(stderr, "-- test_post_build_validation_failure_api\n");
    reset_index_dir();

    CHECK_EQ(run_ikafssnindex(""), 0);

    std::string fp = final_prefix();
    flip_byte(fp + ".kix", /*off=*/96, /*xor_mask=*/0x10);

    Logger silent(Logger::kError);
    CHECK(!validate_volume_final_strict(fp, /*skip_kpx=*/true, silent));
}

// The posting file scratch never survives a completed build, and a stale
// one left by an interrupted run is discarded rather than mistaken for
// resumable postings.
void test_posting_file_scratch_removed() {
    std::fprintf(stderr, "-- test_posting_file_scratch_removed\n");
    reset_index_dir();

    CHECK_EQ(run_ikafssnindex(""), 0);
    CHECK_EQ(count_pf_tmp(), 0);

    // Simulate a crash partway through the postings pass: .kix.pf.tmp
    // holds bytes that belong to no index.
    std::string fp = final_prefix();
    std::string pf = fp + ".kix.pf.tmp";
    write_garbage(pf);
    std::filesystem::remove(fp + ".kix");

    std::string out;
    CHECK_EQ(run_ikafssnindex("", &out), 0);
    CHECK(out.find("=== Writing postings") != std::string::npos);
    CHECK_EQ(count_pf_tmp(), 0);

    Logger silent(Logger::kError);
    CHECK(validate_volume_final_strict(fp, /*skip_kpx=*/true, silent));

    // -force_rebuild 1 sweeps the scratch even when the final files are
    // intact, so it cannot leak into the rebuild.
    write_garbage(pf);
    CHECK_EQ(run_ikafssnindex("-force_rebuild 1"), 0);
    CHECK_EQ(count_pf_tmp(), 0);
    CHECK(validate_volume_final_strict(fp, /*skip_kpx=*/true, silent));
}

void cleanup() {
    std::error_code ec;
    std::filesystem::remove_all(g_index_dir, ec);
}

} // namespace

int main() {
    setup();
    test_resume_full_complete_skip();
    test_resume_one_volume_missing();
    test_resume_one_volume_corrupted();
    test_resume_metadata_tmp_only();
    test_force_rebuild_overrides();
    test_resume_fractional_freq_filter();
    test_resume_khx_only_invalid();
    test_post_build_validation_failure_api();
    test_posting_file_scratch_removed();
    cleanup();
    TEST_SUMMARY();
    return g_fail_count == 0 ? 0 : 1;
}
