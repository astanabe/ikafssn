#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace ikafssn {

// Per-job query file storage for ikafssnhttpd.
//
// Layout: <root>/<ab>/<job_id>.json.zst
//   where <ab> = the first two characters of <job_id> (or padded).
// The file body is always a single zstd frame whose plaintext is the
// JSON request body (the same shape that POST /api/v1/jobs accepts).
// `protocol::serialize(SearchRequest)` is no longer involved in the
// queue persistence path; that binary form is reserved for the wire
// (httpd -> backend ikafssnserver) only.
//
// Two write entry points are provided:
//   - write_compressed_passthrough: store the bytes verbatim.  Used
//     when the client uploaded a Content-Type: application/zstd body;
//     we keep the client's original frame (and its compression level)
//     instead of decompressing and re-compressing.
//   - write_plain: zstd-compress the bytes at `level_` then store.
//     Used when the client uploaded plaintext JSON
//     (Content-Type: application/json or empty).
//
// Writes are atomic (`<path>.tmp` + fdatasync + rename).  On startup,
// any leftover `*.json.zst.tmp` files under <root> are removed by
// `sweep_tmp()` (called from `init()`).
//
// All methods are thread-safe with respect to each other through the
// underlying filesystem (no shared in-process state beyond the immutable
// root_/level_ pair).
class QueryStore {
public:
    QueryStore(std::string root, int level);

    // Create the root directory (mkdir -p) and sweep any leftover
    // *.json.zst.tmp files from a previous unclean shutdown.  Returns
    // false on failure (e.g. permission denied) with `error_msg` set.
    bool init(std::string& error_msg);

    // Persist a pre-compressed zstd frame verbatim (no re-compression).
    // The caller must guarantee `bytes` is a valid zstd frame.
    bool write_compressed_passthrough(const std::string& job_id,
                                      const void* bytes, std::size_t n,
                                      std::string& error_msg);

    // Compress `bytes` with the configured level and persist atomically.
    bool write_plain(const std::string& job_id,
                     const void* bytes, std::size_t n,
                     std::string& error_msg);

    // Read the per-job file and zstd-decompress to `out_json`.  Returns
    // false (with `error_msg` set) if the file is missing, unreadable,
    // or fails to decompress.
    bool read(const std::string& job_id,
              std::string& out_json,
              std::string& error_msg) const;

    // Compute the on-disk path (does not check existence).
    std::string path_for(const std::string& job_id) const;

    // True if the per-job query file exists.
    bool exists(const std::string& job_id) const;

    // Remove the per-job query file.  Missing file is treated as
    // success.  Returns false only on filesystem errors other than
    // "not found".
    bool unlink(const std::string& job_id, std::string& error_msg);

    // Walk the entire root tree and remove every *.json.zst.tmp file
    // (orphans from a crash or hard kill).  Returns the count of files
    // removed; sets error_msg on iteration failure.
    int sweep_tmp(std::string& error_msg);

    // Enumerate every job_id whose <ab>/<job_id>.json.zst file currently
    // exists under <root>.  Used by JobHousekeeper's orphan sweep to
    // find files whose SQLite row was already deleted.
    std::vector<std::string> list_job_ids(std::string& error_msg) const;

    int level() const { return level_; }
    const std::string& root() const { return root_; }

private:
    std::string root_;
    int         level_;

    static std::string shard_for(const std::string& job_id);

    // Atomic writer shared by both write entry points.  `compressed`
    // already holds the bytes that should land verbatim on disk.
    bool write_atomic_(const std::string& job_id,
                       const void* compressed, std::size_t n,
                       std::string& error_msg);
};

} // namespace ikafssn
