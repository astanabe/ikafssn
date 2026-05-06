// Unit test for src/ikafssnhttpd/job_store.cpp.
//
// Covers insert -> fetch -> mark_done round-trip, mark_failed,
// requeue_for_retry (attempts++), delete_expired boundary, and
// requeue_orphans on startup recovery.  Schema-version 3 (file-based
// query+result stores) leaves request and result bodies out of the
// SQLite row entirely; these tests therefore only verify metadata
// transitions and do not assert on body content.

#include "ikafssnhttpd/job_store.hpp"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <thread>
#include <vector>

#include <sqlite3.h>

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

static void test_insert_fetch_done() {
    JobStore s;
    std::string err;
    CHECK(s.open(tmp_db_path(), err));

    int64_t submitted_at = 0;
    bool dup = false;
    CHECK(s.insert_job("job-A", "nt", 3, submitted_at, dup, err));
    CHECK(!dup);
    CHECK(submitted_at > 0);

    JobMeta meta;
    CHECK(s.fetch_one_queued(meta, err));
    CHECK(meta.job_id == "job-A");
    CHECK(meta.status == JobStatus::kRunning);

    // No more queued jobs.
    JobMeta meta2;
    CHECK(!s.fetch_one_queued(meta2, err));

    CHECK(s.mark_done("job-A", err));

    JobMeta got;
    CHECK(s.get_status("job-A", got, err));
    CHECK(got.status == JobStatus::kDone);
    CHECK(got.completed_at > 0);
}

static void test_duplicate_insert() {
    JobStore s;
    std::string err;
    CHECK(s.open(tmp_db_path(), err));

    int64_t at = 0;
    bool dup = false;
    CHECK(s.insert_job("dup", "db", 1, at, dup, err));
    CHECK(!s.insert_job("dup", "db", 1, at, dup, err));
    CHECK(dup);
}

static void test_requeue_retry() {
    JobStore s;
    std::string err;
    CHECK(s.open(tmp_db_path(), err));

    int64_t at = 0;
    bool dup = false;
    CHECK(s.insert_job("retry", "db", 2, at, dup, err));

    JobMeta meta;
    CHECK(s.fetch_one_queued(meta, err));
    CHECK(meta.attempts == 0);

    CHECK(s.requeue_for_retry("retry", "transient", err));

    CHECK(s.fetch_one_queued(meta, err));
    CHECK(meta.attempts == 1);
    CHECK(meta.error_message == "transient");

    CHECK(s.requeue_for_retry("retry", "again", err));
    CHECK(s.fetch_one_queued(meta, err));
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
    CHECK(s.insert_job("old", "db", 1, at, dup, err));
    JobMeta meta;
    CHECK(s.fetch_one_queued(meta, err));
    CHECK(s.mark_done("old", err));

    // Sleep past the retention window.  now_unix() has 1-second
    // granularity so we need >1s of separation to be sure
    // `completed_at < now - retention` holds.
    std::this_thread::sleep_for(std::chrono::milliseconds(2200));

    std::vector<std::string> deleted_ids;
    int64_t n = s.delete_expired(/*retention_seconds=*/1, deleted_ids, err);
    CHECK(n >= 1);
    CHECK(!deleted_ids.empty());
    CHECK(deleted_ids[0] == "old");

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
        CHECK(s.insert_job("orphan", "db", 1, at, dup, err));
        JobMeta meta;
        CHECK(s.fetch_one_queued(meta, err)); // status=running
        // Simulate crash: do not mark_done / mark_failed.
    }
    {
        JobStore s;
        std::string err;
        CHECK(s.open(path, err));
        int64_t n = s.requeue_orphans(err);
        CHECK(n >= 1);

        JobMeta meta;
        CHECK(s.fetch_one_queued(meta, err));
        CHECK(meta.job_id == "orphan");
        CHECK(meta.status == JobStatus::kRunning);
    }
}

// Verifies check_schema() rejects a database that still carries either
// legacy column (result_blob from <v2 or request_blob from v2).
// Synthesises one with raw sqlite3 calls (we cannot use JobStore::open
// since that's exactly what we're testing).
static void test_legacy_schema_rejected_result_blob() {
    auto dir = fs::temp_directory_path() / "ikafssn_job_store_legacy_v1";
    std::error_code ec;
    fs::remove_all(dir, ec);
    fs::create_directories(dir, ec);
    std::string path = (dir / "legacy.db").string();

    sqlite3* h = nullptr;
    int rc = sqlite3_open_v2(path.c_str(), &h,
                             SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE
                                 | SQLITE_OPEN_FULLMUTEX, nullptr);
    CHECK(rc == SQLITE_OK);
    char* eptr = nullptr;
    sqlite3_exec(h,
        "CREATE TABLE jobs ("
        "  job_id TEXT PRIMARY KEY,"
        "  status TEXT NOT NULL,"
        "  result_blob BLOB,"
        "  submitted_at INTEGER NOT NULL"
        ");",
        nullptr, nullptr, &eptr);
    if (eptr) sqlite3_free(eptr);
    sqlite3_exec(h, "PRAGMA user_version=1;", nullptr, nullptr, nullptr);
    sqlite3_close(h);

    JobStore s;
    std::string err;
    CHECK(!s.open(path, err));
    CHECK(err.find("incompatible schema") != std::string::npos);
}

static void test_legacy_schema_rejected_request_blob() {
    auto dir = fs::temp_directory_path() / "ikafssn_job_store_legacy_v2";
    std::error_code ec;
    fs::remove_all(dir, ec);
    fs::create_directories(dir, ec);
    std::string path = (dir / "legacy.db").string();

    sqlite3* h = nullptr;
    int rc = sqlite3_open_v2(path.c_str(), &h,
                             SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE
                                 | SQLITE_OPEN_FULLMUTEX, nullptr);
    CHECK(rc == SQLITE_OK);
    char* eptr = nullptr;
    sqlite3_exec(h,
        "CREATE TABLE jobs ("
        "  job_id TEXT PRIMARY KEY,"
        "  status TEXT NOT NULL,"
        "  request_blob BLOB NOT NULL,"
        "  submitted_at INTEGER NOT NULL"
        ");",
        nullptr, nullptr, &eptr);
    if (eptr) sqlite3_free(eptr);
    // user_version=2 corresponds to the previous (request_blob in row)
    // schema; with user_version=2 + has_request_blob the new code must
    // still reject it.
    sqlite3_exec(h, "PRAGMA user_version=2;", nullptr, nullptr, nullptr);
    sqlite3_close(h);

    JobStore s;
    std::string err;
    CHECK(!s.open(path, err));
    CHECK(err.find("incompatible schema") != std::string::npos);
    CHECK(err.find("request_blob") != std::string::npos);
}

int main() {
    test_insert_fetch_done();
    test_duplicate_insert();
    test_requeue_retry();
    test_delete_expired();
    test_requeue_orphans();
    test_legacy_schema_rejected_result_blob();
    test_legacy_schema_rejected_request_blob();
    TEST_SUMMARY();
    return g_fail_count == 0 ? 0 : 1;
}
