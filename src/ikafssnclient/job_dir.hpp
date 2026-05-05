#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "io/result_writer.hpp"

namespace ikafssn {

// Per-job persistent state at .ikafssnclient/<group_id>/<job_id>.json.
struct JobMeta {
    std::string job_id;
    std::string group_id;
    int32_t     n_seqs = 0;
    std::string fasta_file;        // canonical path to query FASTA
    std::pair<int32_t, int32_t> seq_index_range = {0, 0};

    std::string status;            // "submitted","queued","running","done","failed","timeout"
    int32_t     attempts = 0;
    int64_t     submitted_at = 0;
    int64_t     completed_at = 0;
    int64_t     last_polled_at = 0;
    std::string error_message;
    std::string fail_reason;
};

// Per-group metadata at .ikafssnclient/<group_id>/group.json.
struct GroupMeta {
    std::string group_id;
    int64_t     submitted_at = 0;
    std::string httpd_url;
    std::string db;
    std::string query_file_path_abs;
    std::string query_file_sha256;
    int32_t     max_nseq_per_req = 0;
    uint8_t     k = 0;
    uint8_t     mode = 0;
    uint8_t     t = 0;
    uint8_t     template_type = 0;
    std::string output_format;
    std::string output_path;
    int32_t     compression_level = -1;  // -1 = codec default
    std::vector<std::string> job_ids;

    int32_t cnt_queued = 0;
    int32_t cnt_running = 0;
    int32_t cnt_done = 0;
    int32_t cnt_failed = 0;
};

// Top-level helpers operating on the `.ikafssnclient/` directory tree.

std::string default_jobs_root();
bool        ensure_jobs_root(std::string& error_msg);

// Persist GroupMeta to <root>/<group_id>/group.json (atomic tmp+rename).
bool write_group_meta(const std::string& root, const GroupMeta& meta,
                      std::string& error_msg);
bool read_group_meta(const std::string& root, const std::string& group_id,
                     GroupMeta& out, std::string& error_msg);

// Persist JobMeta to <root>/<group_id>/<job_id>.json (atomic tmp+rename).
bool write_job_meta(const std::string& root, const JobMeta& meta,
                    std::string& error_msg);
bool read_job_meta(const std::string& root, const std::string& group_id,
                   const std::string& job_id, JobMeta& out,
                   std::string& error_msg);

// List all group_ids under <root>, ordered by GroupMeta.submitted_at.
std::vector<std::string> list_groups(const std::string& root);

// List all job_ids inside a single group.
std::vector<std::string> list_jobs(const std::string& root,
                                   const std::string& group_id);

// Given an arbitrary id supplied to -resume / -job_detail, work out
// whether it names a group or a job (and if a job, which group).  On
// success returns true; sets `is_group=true` and `group_id=id` for a
// group, or `is_group=false` + group_id of the parent for a job.
bool resolve_id(const std::string& root, const std::string& id,
                bool& is_group, std::string& group_id, std::string& job_id_out);

// Defline cache (`<job_id>.deflines.zst`).  One LF-terminated line per
// query in seq_index_range.  Used so failed_writer can emit
// *FAILED:<reason> rows even when the original FASTA is no longer on disk.
bool write_job_deflines(const std::string& root,
                        const std::string& group_id,
                        const std::string& job_id,
                        const std::vector<std::string>& deflines,
                        std::string& error_msg);
bool read_job_deflines(const std::string& root,
                       const std::string& group_id,
                       const std::string& job_id,
                       std::vector<std::string>& out,
                       std::string& error_msg);

// Result blob cache (`<job_id>.result.bin.zst`).  Stores the raw
// serialize(SearchResponse) payload returned by GET /jobs/<id>/result
// before it has been merged into the user's output file.
bool write_job_result(const std::string& root,
                      const std::string& group_id,
                      const std::string& job_id,
                      const std::vector<uint8_t>& blob,
                      std::string& error_msg);
bool read_job_result(const std::string& root,
                     const std::string& group_id,
                     const std::string& job_id,
                     std::vector<uint8_t>& out,
                     std::string& error_msg);
bool delete_job_result(const std::string& root,
                       const std::string& group_id,
                       const std::string& job_id,
                       std::string& error_msg);

// Generate a fresh UUIDv4 string using a process-local RNG.
std::string make_uuidv4();

} // namespace ikafssn
