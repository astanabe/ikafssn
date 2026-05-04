// Unit test for src/ikafssnhttpd/job_store.cpp.
//
// Covers insert -> fetch -> mark_done round-trip, mark_failed,
// requeue_for_retry (attempts++), delete_expired boundary, and
// requeue_orphans on startup recovery.

#include "ikafssnhttpd/job_store.hpp"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <thread>
#include <vector>

#include "test_util.hpp"

namespace fs = std::filesystem;

using ikafssn::JobMeta;
using ikafssn::JobStatus;
using ikafssn::JobStore;

static std::string tmp_db_path() {
    auto dir = fs::temp_directory_path() / "ikafssn_job_store_test";
    std::error_code ec;
    fs::remove_all(dir, ec);
    fs::create_directories(dir, ec);
    return (dir / "jobs.db").string();
}

static std::vector<uint8_t> blob(const std::string& s) {
    return std::vector<uint8_t>(s.begin(), s.end());
}

static void test_insert_fetch_done() {
    JobStore s;
    std::string err;
    CHECK(s.open(tmp_db_path(), err));

    int64_t submitted_at = 0;
    bool dup = false;
    CHECK(s.insert_job("job-A", "nt", 3, blob("REQA"),
                       submitted_at, dup, err));
    CHECK(!dup);
    CHECK(submitted_at > 0);

    JobMeta meta;
    std::vector<uint8_t> req;
    CHECK(s.fetch_one_queued(meta, req, err));
    CHECK(meta.job_id == "job-A");
    CHECK(meta.status == JobStatus::kRunning);
    CHECK(req.size() == 4);

    // No more queued jobs.
    JobMeta meta2;
    std::vector<uint8_t> req2;
    CHECK(!s.fetch_one_queued(meta2, req2, err));

    CHECK(s.mark_done("job-A", blob("RESPA"), err));

    JobMeta got;
    CHECK(s.get_status("job-A", got, err));
    CHECK(got.status == JobStatus::kDone);

    bool nf = false, ws = false;
    std::vector<uint8_t> result;
    CHECK(s.get_result("job-A", result, nf, ws, err));
    CHECK(!nf);
    CHECK(!ws);
    CHECK(result.size() == 5);
}

static void test_duplicate_insert() {
    JobStore s;
    std::string err;
    CHECK(s.open(tmp_db_path(), err));

    int64_t at = 0;
    bool dup = false;
    CHECK(s.insert_job("dup", "db", 1, blob("X"), at, dup, err));
    CHECK(!s.insert_job("dup", "db", 1, blob("X"), at, dup, err));
    CHECK(dup);
}

static void test_requeue_retry() {
    JobStore s;
    std::string err;
    CHECK(s.open(tmp_db_path(), err));

    int64_t at = 0;
    bool dup = false;
    CHECK(s.insert_job("retry", "db", 2, blob("R"), at, dup, err));

    JobMeta meta;
    std::vector<uint8_t> req;
    CHECK(s.fetch_one_queued(meta, req, err));
    CHECK(meta.attempts == 0);

    CHECK(s.requeue_for_retry("retry", "transient", err));

    CHECK(s.fetch_one_queued(meta, req, err));
    CHECK(meta.attempts == 1);
    CHECK(meta.error_message == "transient");

    CHECK(s.requeue_for_retry("retry", "again", err));
    CHECK(s.fetch_one_queued(meta, req, err));
    CHECK(meta.attempts == 2);

    CHECK(s.mark_failed("retry", "final", "backend_unreachable: x", err));
    JobMeta got;
    CHECK(s.get_status("retry", got, err));
    CHECK(got.status == JobStatus::kFailed);
    CHECK(got.fail_reason == "backend_unreachable: x");
}

static void test_delete_expired() {
    JobStore s;
    std::string err;
    CHECK(s.open(tmp_db_path(), err));

    int64_t at = 0;
    bool dup = false;
    CHECK(s.insert_job("old", "db", 1, blob("o"), at, dup, err));
    JobMeta meta; std::vector<uint8_t> req;
    CHECK(s.fetch_one_queued(meta, req, err));
    CHECK(s.mark_done("old", blob("r"), err));

    // Sleep past the retention window.  now_unix() has 1-second
    // granularity so we need >1s of separation to be sure
    // `completed_at < now - retention` holds.
    std::this_thread::sleep_for(std::chrono::milliseconds(2200));

    int64_t n = s.delete_expired(/*retention_seconds=*/1, err);
    CHECK(n >= 1);

    JobMeta got;
    CHECK(!s.get_status("old", got, err));
}

static void test_requeue_orphans() {
    auto path = tmp_db_path();
    {
        JobStore s;
        std::string err;
        CHECK(s.open(path, err));
        int64_t at = 0; bool dup = false;
        CHECK(s.insert_job("orphan", "db", 1, blob("o"), at, dup, err));
        JobMeta meta; std::vector<uint8_t> req;
        CHECK(s.fetch_one_queued(meta, req, err)); // status=running
        // Simulate crash: do not mark_done / mark_failed.
    }
    {
        JobStore s;
        std::string err;
        CHECK(s.open(path, err));
        int64_t n = s.requeue_orphans(err);
        CHECK(n >= 1);

        JobMeta meta; std::vector<uint8_t> req;
        CHECK(s.fetch_one_queued(meta, req, err));
        CHECK(meta.job_id == "orphan");
        CHECK(meta.status == JobStatus::kRunning);
    }
}

int main() {
    test_insert_fetch_done();
    test_duplicate_insert();
    test_requeue_retry();
    test_delete_expired();
    test_requeue_orphans();
    TEST_SUMMARY();
    return g_fail_count == 0 ? 0 : 1;
}
