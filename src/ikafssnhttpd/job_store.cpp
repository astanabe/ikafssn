#include "ikafssnhttpd/job_store.hpp"

#include <chrono>
#include <cstring>
#include <filesystem>

#include <sqlite3.h>

namespace ikafssn {

namespace {

int64_t now_unix() {
    using namespace std::chrono;
    return duration_cast<seconds>(system_clock::now().time_since_epoch()).count();
}

} // namespace

const char* job_status_str(JobStatus s) {
    switch (s) {
        case JobStatus::kQueued:  return "queued";
        case JobStatus::kRunning: return "running";
        case JobStatus::kDone:    return "done";
        case JobStatus::kFailed:  return "failed";
    }
    return "queued";
}

bool job_status_parse(const std::string& s, JobStatus& out) {
    if (s == "queued")  { out = JobStatus::kQueued;  return true; }
    if (s == "running") { out = JobStatus::kRunning; return true; }
    if (s == "done")    { out = JobStatus::kDone;    return true; }
    if (s == "failed")  { out = JobStatus::kFailed;  return true; }
    return false;
}

JobStore::~JobStore() { close(); }

bool JobStore::exec_(const char* sql, std::string& error_msg) {
    char* err = nullptr;
    int rc = sqlite3_exec(db_, sql, nullptr, nullptr, &err);
    if (rc != SQLITE_OK) {
        error_msg = std::string("sqlite exec failed: ")
                  + (err ? err : "unknown");
        if (err) sqlite3_free(err);
        return false;
    }
    return true;
}

bool JobStore::check_schema(std::string& error_msg) {
    // user_version probe.
    sqlite3_stmt* stmt = nullptr;
    int rc = sqlite3_prepare_v2(db_, "PRAGMA user_version;", -1,
                                &stmt, nullptr);
    if (rc != SQLITE_OK) {
        error_msg = std::string("prepare PRAGMA user_version: ")
                  + sqlite3_errmsg(db_);
        return false;
    }
    int user_version = 0;
    if (sqlite3_step(stmt) == SQLITE_ROW) {
        user_version = sqlite3_column_int(stmt, 0);
    }
    sqlite3_finalize(stmt);

    // Probe whether the jobs table exists at all and whether it carries
    // the legacy `result_blob` column.
    bool has_jobs       = false;
    bool has_result_blob = false;
    rc = sqlite3_prepare_v2(db_, "PRAGMA table_info(jobs);", -1,
                            &stmt, nullptr);
    if (rc != SQLITE_OK) {
        error_msg = std::string("prepare PRAGMA table_info: ")
                  + sqlite3_errmsg(db_);
        return false;
    }
    while (sqlite3_step(stmt) == SQLITE_ROW) {
        has_jobs = true;
        const unsigned char* name = sqlite3_column_text(stmt, 1);
        if (name && std::strcmp(reinterpret_cast<const char*>(name),
                                "result_blob") == 0) {
            has_result_blob = true;
        }
    }
    sqlite3_finalize(stmt);

    if (!has_jobs) {
        // Brand-new DB; open() will create the schema and stamp
        // user_version=2.
        return true;
    }
    if (has_result_blob || user_version < 2) {
        error_msg = "jobs.db has incompatible schema (user_version="
                  + std::to_string(user_version)
                  + (has_result_blob ? ", legacy result_blob column present"
                                     : "")
                  + "); please delete jobs.db and restart so the new"
                    " file-based result store schema can be initialised";
        return false;
    }
    return true;
}

bool JobStore::open(const std::string& path, std::string& error_msg) {
    std::lock_guard<std::mutex> lock(mu_);
    if (db_) {
        error_msg = "JobStore already open";
        return false;
    }

    // Ensure parent directory exists.
    std::filesystem::path p(path);
    if (p.has_parent_path()) {
        std::error_code ec;
        std::filesystem::create_directories(p.parent_path(), ec);
        // ignore EEXIST; error reported by sqlite3_open_v2 if it really fails.
    }

    int rc = sqlite3_open_v2(path.c_str(), &db_,
                             SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE
                                 | SQLITE_OPEN_FULLMUTEX,
                             nullptr);
    if (rc != SQLITE_OK) {
        error_msg = "sqlite3_open_v2 failed for " + path + ": ";
        error_msg += db_ ? sqlite3_errmsg(db_) : "unknown";
        if (db_) { sqlite3_close(db_); db_ = nullptr; }
        return false;
    }

    // PRAGMAs as required by the design doc.
    if (!exec_("PRAGMA journal_mode=WAL;",     error_msg)) return false;
    if (!exec_("PRAGMA synchronous=NORMAL;",   error_msg)) return false;
    if (!exec_("PRAGMA busy_timeout=5000;",    error_msg)) return false;

    if (!check_schema(error_msg)) {
        sqlite3_close(db_);
        db_ = nullptr;
        return false;
    }

    if (!exec_(
        "CREATE TABLE IF NOT EXISTS jobs ("
        " job_id        TEXT     PRIMARY KEY,"
        " status        TEXT     NOT NULL CHECK ("
        "   status IN ('queued','running','done','failed')),"
        " request_blob  BLOB     NOT NULL,"
        " error_message TEXT,"
        " fail_reason   TEXT,"
        " attempts      INTEGER  NOT NULL DEFAULT 0,"
        " submitted_at  INTEGER  NOT NULL,"
        " started_at    INTEGER,"
        " completed_at  INTEGER,"
        " db            TEXT,"
        " n_seqs        INTEGER"
        ");",
        error_msg)) return false;
    if (!exec_("CREATE INDEX IF NOT EXISTS idx_jobs_status_submitted"
               " ON jobs(status, submitted_at);", error_msg)) return false;
    if (!exec_("CREATE INDEX IF NOT EXISTS idx_jobs_completed_at"
               " ON jobs(completed_at);", error_msg)) return false;

    // Stamp user_version=2 (file-based result store).  Do this last so
    // a partially-initialised DB still fails check_schema on the next
    // startup.
    if (!exec_("PRAGMA user_version=2;", error_msg)) return false;

    return true;
}

void JobStore::close() {
    std::lock_guard<std::mutex> lock(mu_);
    if (db_) {
        sqlite3_close(db_);
        db_ = nullptr;
    }
}

bool JobStore::insert_job(const std::string& job_id,
                          const std::string& db,
                          int32_t n_seqs,
                          const std::vector<uint8_t>& request_blob,
                          int64_t& submitted_at,
                          bool& duplicate,
                          std::string& error_msg) {
    std::lock_guard<std::mutex> lock(mu_);
    duplicate = false;
    submitted_at = now_unix();

    sqlite3_stmt* stmt = nullptr;
    const char* sql =
        "INSERT INTO jobs(job_id, status, request_blob, attempts,"
        "                 submitted_at, db, n_seqs)"
        " VALUES(?, 'queued', ?, 0, ?, ?, ?);";
    int rc = sqlite3_prepare_v2(db_, sql, -1, &stmt, nullptr);
    if (rc != SQLITE_OK) {
        error_msg = std::string("prepare insert: ") + sqlite3_errmsg(db_);
        return false;
    }
    sqlite3_bind_text(stmt, 1, job_id.data(),
                      static_cast<int>(job_id.size()), SQLITE_TRANSIENT);
    sqlite3_bind_blob(stmt, 2, request_blob.data(),
                      static_cast<int>(request_blob.size()), SQLITE_TRANSIENT);
    sqlite3_bind_int64(stmt, 3, submitted_at);
    sqlite3_bind_text(stmt, 4, db.data(),
                      static_cast<int>(db.size()), SQLITE_TRANSIENT);
    sqlite3_bind_int(stmt, 5, n_seqs);
    rc = sqlite3_step(stmt);
    sqlite3_finalize(stmt);

    if (rc == SQLITE_CONSTRAINT) {
        duplicate = true;
        return false;
    }
    if (rc != SQLITE_DONE) {
        error_msg = std::string("insert step: ") + sqlite3_errmsg(db_);
        return false;
    }
    return true;
}

bool JobStore::fetch_one_queued(JobMeta& out_meta,
                                std::vector<uint8_t>& out_request_blob,
                                std::string& error_msg) {
    std::lock_guard<std::mutex> lock(mu_);

    int64_t now = now_unix();
    sqlite3_stmt* stmt = nullptr;
    const char* sql =
        "UPDATE jobs SET status='running', started_at=? "
        " WHERE job_id = ("
        "   SELECT job_id FROM jobs"
        "    WHERE status='queued'"
        "    ORDER BY submitted_at ASC LIMIT 1)"
        " RETURNING job_id, request_blob, attempts, submitted_at,"
        "          started_at, db, n_seqs, error_message;";
    int rc = sqlite3_prepare_v2(db_, sql, -1, &stmt, nullptr);
    if (rc != SQLITE_OK) {
        error_msg = std::string("prepare fetch: ") + sqlite3_errmsg(db_);
        return false;
    }
    sqlite3_bind_int64(stmt, 1, now);

    rc = sqlite3_step(stmt);
    if (rc == SQLITE_DONE) {
        sqlite3_finalize(stmt);
        return false;
    }
    if (rc != SQLITE_ROW) {
        error_msg = std::string("fetch step: ") + sqlite3_errmsg(db_);
        sqlite3_finalize(stmt);
        return false;
    }

    out_meta = JobMeta{};
    out_meta.status = JobStatus::kRunning;

    const unsigned char* jid = sqlite3_column_text(stmt, 0);
    if (jid) out_meta.job_id.assign(reinterpret_cast<const char*>(jid));

    const void* blob = sqlite3_column_blob(stmt, 1);
    int blob_n = sqlite3_column_bytes(stmt, 1);
    out_request_blob.assign(reinterpret_cast<const uint8_t*>(blob),
                            reinterpret_cast<const uint8_t*>(blob) + blob_n);

    out_meta.attempts     = sqlite3_column_int(stmt, 2);
    out_meta.submitted_at = sqlite3_column_int64(stmt, 3);
    out_meta.started_at   = sqlite3_column_int64(stmt, 4);
    const unsigned char* db_text = sqlite3_column_text(stmt, 5);
    if (db_text) out_meta.db.assign(reinterpret_cast<const char*>(db_text));
    out_meta.n_seqs       = sqlite3_column_int(stmt, 6);
    const unsigned char* em_text = sqlite3_column_text(stmt, 7);
    if (em_text) out_meta.error_message.assign(
        reinterpret_cast<const char*>(em_text));

    sqlite3_finalize(stmt);
    return true;
}

bool JobStore::mark_done(const std::string& job_id, std::string& error_msg) {
    std::lock_guard<std::mutex> lock(mu_);
    int64_t now = now_unix();

    sqlite3_stmt* stmt = nullptr;
    const char* sql =
        "UPDATE jobs SET status='done', completed_at=?"
        " WHERE job_id=? AND status='running';";
    int rc = sqlite3_prepare_v2(db_, sql, -1, &stmt, nullptr);
    if (rc != SQLITE_OK) {
        error_msg = std::string("prepare mark_done: ") + sqlite3_errmsg(db_);
        return false;
    }
    sqlite3_bind_int64(stmt, 1, now);
    sqlite3_bind_text(stmt, 2, job_id.data(),
                      static_cast<int>(job_id.size()), SQLITE_TRANSIENT);
    rc = sqlite3_step(stmt);
    sqlite3_finalize(stmt);
    if (rc != SQLITE_DONE) {
        error_msg = std::string("mark_done step: ") + sqlite3_errmsg(db_);
        return false;
    }
    return true;
}

bool JobStore::mark_failed(const std::string& job_id,
                           const std::string& error_message,
                           const std::string& fail_reason,
                           std::string& error_msg) {
    std::lock_guard<std::mutex> lock(mu_);
    int64_t now = now_unix();

    sqlite3_stmt* stmt = nullptr;
    const char* sql =
        "UPDATE jobs SET status='failed', error_message=?, fail_reason=?,"
        "                completed_at=?"
        " WHERE job_id=?;";
    int rc = sqlite3_prepare_v2(db_, sql, -1, &stmt, nullptr);
    if (rc != SQLITE_OK) {
        error_msg = std::string("prepare mark_failed: ") + sqlite3_errmsg(db_);
        return false;
    }
    sqlite3_bind_text(stmt, 1, error_message.data(),
                      static_cast<int>(error_message.size()), SQLITE_TRANSIENT);
    sqlite3_bind_text(stmt, 2, fail_reason.data(),
                      static_cast<int>(fail_reason.size()), SQLITE_TRANSIENT);
    sqlite3_bind_int64(stmt, 3, now);
    sqlite3_bind_text(stmt, 4, job_id.data(),
                      static_cast<int>(job_id.size()), SQLITE_TRANSIENT);
    rc = sqlite3_step(stmt);
    sqlite3_finalize(stmt);
    if (rc != SQLITE_DONE) {
        error_msg = std::string("mark_failed step: ") + sqlite3_errmsg(db_);
        return false;
    }
    return true;
}

bool JobStore::requeue_for_retry(const std::string& job_id,
                                 const std::string& error_message,
                                 std::string& error_msg) {
    std::lock_guard<std::mutex> lock(mu_);

    sqlite3_stmt* stmt = nullptr;
    const char* sql =
        "UPDATE jobs SET status='queued',"
        "                attempts = attempts + 1,"
        "                error_message = ?"
        " WHERE job_id=? AND status='running';";
    int rc = sqlite3_prepare_v2(db_, sql, -1, &stmt, nullptr);
    if (rc != SQLITE_OK) {
        error_msg = std::string("prepare requeue: ") + sqlite3_errmsg(db_);
        return false;
    }
    sqlite3_bind_text(stmt, 1, error_message.data(),
                      static_cast<int>(error_message.size()), SQLITE_TRANSIENT);
    sqlite3_bind_text(stmt, 2, job_id.data(),
                      static_cast<int>(job_id.size()), SQLITE_TRANSIENT);
    rc = sqlite3_step(stmt);
    sqlite3_finalize(stmt);
    if (rc != SQLITE_DONE) {
        error_msg = std::string("requeue step: ") + sqlite3_errmsg(db_);
        return false;
    }
    return true;
}

bool JobStore::get_status(const std::string& job_id, JobMeta& out,
                          std::string& error_msg) {
    std::lock_guard<std::mutex> lock(mu_);

    sqlite3_stmt* stmt = nullptr;
    const char* sql =
        "SELECT status, error_message, fail_reason, attempts,"
        "       submitted_at, started_at, completed_at, db, n_seqs"
        "  FROM jobs WHERE job_id=?;";
    int rc = sqlite3_prepare_v2(db_, sql, -1, &stmt, nullptr);
    if (rc != SQLITE_OK) {
        error_msg = std::string("prepare get_status: ") + sqlite3_errmsg(db_);
        return false;
    }
    sqlite3_bind_text(stmt, 1, job_id.data(),
                      static_cast<int>(job_id.size()), SQLITE_TRANSIENT);

    rc = sqlite3_step(stmt);
    if (rc == SQLITE_DONE) {
        sqlite3_finalize(stmt);
        return false;
    }
    if (rc != SQLITE_ROW) {
        error_msg = std::string("get_status step: ") + sqlite3_errmsg(db_);
        sqlite3_finalize(stmt);
        return false;
    }

    out = JobMeta{};
    out.job_id = job_id;
    const unsigned char* st = sqlite3_column_text(stmt, 0);
    if (st) {
        std::string s(reinterpret_cast<const char*>(st));
        job_status_parse(s, out.status);
    }
    auto col_text = [&](int i, std::string& dst) {
        const unsigned char* t = sqlite3_column_text(stmt, i);
        if (t) dst.assign(reinterpret_cast<const char*>(t));
    };
    col_text(1, out.error_message);
    col_text(2, out.fail_reason);
    out.attempts     = sqlite3_column_int(stmt, 3);
    out.submitted_at = sqlite3_column_int64(stmt, 4);
    out.started_at   = sqlite3_column_int64(stmt, 5);
    out.completed_at = sqlite3_column_int64(stmt, 6);
    col_text(7, out.db);
    out.n_seqs       = sqlite3_column_int(stmt, 8);

    sqlite3_finalize(stmt);
    return true;
}

int64_t JobStore::delete_expired(int64_t retention_seconds,
                                 std::vector<std::string>& out_deleted_ids,
                                 std::string& error_msg) {
    std::lock_guard<std::mutex> lock(mu_);
    int64_t cutoff = now_unix() - retention_seconds;
    int64_t deleted = 0;

    if (!exec_("BEGIN IMMEDIATE;", error_msg)) return 0;

    {
        sqlite3_stmt* stmt = nullptr;
        const char* sql =
            "SELECT job_id FROM jobs"
            " WHERE status IN ('done','failed')"
            "   AND completed_at IS NOT NULL"
            "   AND completed_at < ?;";
        int rc = sqlite3_prepare_v2(db_, sql, -1, &stmt, nullptr);
        if (rc != SQLITE_OK) {
            error_msg = std::string("prepare delete_expired select: ")
                      + sqlite3_errmsg(db_);
            std::string ignored;
            exec_("ROLLBACK;", ignored);
            return 0;
        }
        sqlite3_bind_int64(stmt, 1, cutoff);
        while ((rc = sqlite3_step(stmt)) == SQLITE_ROW) {
            const unsigned char* jid = sqlite3_column_text(stmt, 0);
            if (jid) out_deleted_ids.emplace_back(
                reinterpret_cast<const char*>(jid));
        }
        sqlite3_finalize(stmt);
        if (rc != SQLITE_DONE) {
            error_msg = std::string("delete_expired select step: ")
                      + sqlite3_errmsg(db_);
            std::string ignored;
            exec_("ROLLBACK;", ignored);
            return 0;
        }
    }

    {
        sqlite3_stmt* stmt = nullptr;
        const char* sql =
            "DELETE FROM jobs"
            " WHERE status IN ('done','failed')"
            "   AND completed_at IS NOT NULL"
            "   AND completed_at < ?;";
        int rc = sqlite3_prepare_v2(db_, sql, -1, &stmt, nullptr);
        if (rc != SQLITE_OK) {
            error_msg = std::string("prepare delete_expired delete: ")
                      + sqlite3_errmsg(db_);
            std::string ignored;
            exec_("ROLLBACK;", ignored);
            return 0;
        }
        sqlite3_bind_int64(stmt, 1, cutoff);
        rc = sqlite3_step(stmt);
        sqlite3_finalize(stmt);
        if (rc != SQLITE_DONE) {
            error_msg = std::string("delete_expired delete step: ")
                      + sqlite3_errmsg(db_);
            std::string ignored;
            exec_("ROLLBACK;", ignored);
            return 0;
        }
        deleted = sqlite3_changes(db_);
    }

    if (!exec_("COMMIT;", error_msg)) {
        std::string ignored;
        exec_("ROLLBACK;", ignored);
        return 0;
    }
    return deleted;
}

int64_t JobStore::requeue_orphans(std::string& error_msg) {
    std::lock_guard<std::mutex> lock(mu_);
    sqlite3_stmt* stmt = nullptr;
    const char* sql =
        "UPDATE jobs SET status='queued', started_at=NULL"
        " WHERE status='running';";
    int rc = sqlite3_prepare_v2(db_, sql, -1, &stmt, nullptr);
    if (rc != SQLITE_OK) {
        error_msg = std::string("prepare requeue_orphans: ") + sqlite3_errmsg(db_);
        return 0;
    }
    rc = sqlite3_step(stmt);
    sqlite3_finalize(stmt);
    if (rc != SQLITE_DONE) {
        error_msg = std::string("requeue_orphans step: ") + sqlite3_errmsg(db_);
        return 0;
    }
    return sqlite3_changes(db_);
}

bool JobStore::prepare_or_load_(const char*, const char*, sqlite3_stmt**,
                                std::string&) {
    // Reserved for a future prepared-statement cache.  All current call
    // sites prepare and finalize per call to keep the implementation
    // simple under the global mutex.
    return false;
}

} // namespace ikafssn
