#pragma once

#include <cstdint>
#include <mutex>
#include <optional>
#include <string>
#include <vector>

struct sqlite3;
struct sqlite3_stmt;

namespace ikafssn {

// Job lifecycle states.  See doc/ikafssn.en.md async REST section.
enum class JobStatus : uint8_t {
    kQueued  = 0,
    kRunning = 1,
    kDone    = 2,
    kFailed  = 3,
};

const char* job_status_str(JobStatus s);
bool        job_status_parse(const std::string& s, JobStatus& out);

// In-memory snapshot of a SQLite row (without the heavy blobs).  Used by
// JobWorker to drive route_search and by HttpController to answer
// GET /api/v1/jobs/<id>.
struct JobMeta {
    std::string job_id;
    JobStatus   status = JobStatus::kQueued;
    std::string error_message;
    std::string fail_reason;
    int32_t     attempts = 0;
    int64_t     submitted_at = 0;
    int64_t     started_at   = 0;
    int64_t     completed_at = 0;
    std::string db;
    int32_t     n_seqs = 0;
};

// SQLite-backed persistent job store for ikafssnhttpd's async REST stack.
//
// All methods serialise through a single internal mutex + a single sqlite3*
// handle; callers do not need to coordinate.  PRAGMAs (WAL / NORMAL /
// busy_timeout=5000) are applied in `open()`.
//
// `fetch_one_queued` uses SQLite 3.35+ `UPDATE ... RETURNING` to pop the
// oldest queued job atomically (Ubuntu 22.04 ships 3.37, 24.04 ships 3.45).
//
// Schema version 3 (current): the SQLite row carries only metadata.
// Both the result body (in ResultStore) and the request body (in
// QueryStore) live as per-job files on disk.  An older schema that
// still has a `result_blob` or `request_blob` column, or has
// user_version<3, is rejected by `check_schema()` so the operator must
// delete jobs.db (and empty the QueryStore directory) before starting
// the new server.
class JobStore {
public:
    JobStore() = default;
    ~JobStore();

    JobStore(const JobStore&) = delete;
    JobStore& operator=(const JobStore&) = delete;

    // Open the SQLite DB and apply PRAGMAs / DDL.  Creates parent dir and
    // file if needed.  Returns false (with `error_msg` set) on failure.
    // Existing databases written by the old schema (containing
    // `result_blob` or `request_blob`, or with `user_version < 3`) are
    // rejected here.
    bool open(const std::string& path, std::string& error_msg);

    // Close the underlying handle.  Idempotent.
    void close();

    // Insert a new job.  The query body is persisted by the caller into
    // the QueryStore before this call (and unlinked on insert failure).
    // Returns false if `job_id` already exists (caller should map that
    // to HTTP 409).  On success, fills `submitted_at` with the row's
    // server-side timestamp.
    bool insert_job(const std::string& job_id,
                    const std::string& db,
                    int32_t n_seqs,
                    int64_t& submitted_at,
                    bool& duplicate,
                    std::string& error_msg);

    // Atomic-pop one queued job (oldest by submitted_at) and flip its
    // status to 'running'.  Sets `out_meta`.  The query body is read
    // from the QueryStore by the caller using `out_meta.job_id`.
    // Returns false (with empty error_msg) when there is no queued job.
    bool fetch_one_queued(JobMeta& out_meta,
                          std::string& error_msg);

    // Mark a job as completed.  The result file is the caller's
    // responsibility to write before calling this — see ResultStore.
    bool mark_done(const std::string& job_id, std::string& error_msg);

    // Mark a job as terminally failed (no further retries).  Records
    // `error_message` and `fail_reason`.
    bool mark_failed(const std::string& job_id,
                     const std::string& error_message,
                     const std::string& fail_reason,
                     std::string& error_msg);

    // Return the job to the queue for another attempt, incrementing
    // `attempts`.  Used when a retryable backend error is seen and at
    // least one backend is still healthy.
    bool requeue_for_retry(const std::string& job_id,
                           const std::string& error_message,
                           std::string& error_msg);

    // Read the meta row.  Returns false if job_id not found.
    bool get_status(const std::string& job_id, JobMeta& out,
                    std::string& error_msg);

    // Delete done/failed rows whose `completed_at` is older than
    // `now_unix - retention_seconds`.  Returns number of rows deleted.
    // The deleted ids are appended to `out_deleted_ids` so the caller
    // (housekeeper) can also remove the matching ResultStore files.
    // SELECT and DELETE happen inside a single transaction so a row can
    // never be reported deleted but still selectable by another thread.
    int64_t delete_expired(int64_t retention_seconds,
                           std::vector<std::string>& out_deleted_ids,
                           std::string& error_msg);

    // On startup, flip any leftover `running` rows back to `queued`.
    // Called once before the worker pool starts taking jobs.
    int64_t requeue_orphans(std::string& error_msg);

    // Inspect `PRAGMA user_version` and `PRAGMA table_info(jobs)` and
    // refuse to run against a schema that predates the file-based
    // query+result stores (user_version<3, or `result_blob` /
    // `request_blob` column present).  Called inside `open()` before
    // any DDL is applied.
    bool check_schema(std::string& error_msg);

private:
    sqlite3*               db_ = nullptr;
    mutable std::mutex     mu_;

    bool exec_(const char* sql, std::string& error_msg);
};

} // namespace ikafssn
