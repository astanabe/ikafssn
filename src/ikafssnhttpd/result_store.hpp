#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace ikafssn {

// Per-job result file storage for ikafssnhttpd.
//
// Layout: <root>/<ab>/<job_id>.bin.zst
//   where <ab> = the first two characters of <job_id> (or padded).
// The file body is a single zstd frame containing the raw
// serialize(SearchResponse) bytes; no additional framing or magic.
//
// Writes are atomic (`<path>.tmp` + fdatasync + rename).  On startup,
// any leftover `*.bin.zst.tmp` files under <root> are removed by
// `sweep_tmp()` (called from `init()`).
//
// All methods are thread-safe with respect to each other through the
// underlying filesystem (no shared in-process state beyond the immutable
// root_/level_ pair).
class ResultStore {
public:
    ResultStore(std::string root, int level);

    // Create the root directory (mkdir -p) and sweep any leftover
    // *.bin.zst.tmp files from a previous unclean shutdown.  Returns
    // false on failure (e.g. permission denied) with `error_msg` set.
    bool init(std::string& error_msg);

    // Write a per-job result file: zstd-compress the bytes at the
    // configured level and persist atomically.  Returns false on
    // failure.  Creates the per-job shard directory if needed.
    bool write(const std::string& job_id,
               const void* bytes, std::size_t n,
               std::string& error_msg);

    // Compute the on-disk path (does not check existence).
    std::string path_for(const std::string& job_id) const;

    // True if the per-job result file exists.
    bool exists(const std::string& job_id) const;

    // Remove the per-job result file.  Missing file is treated as
    // success.  Returns false only on filesystem errors other than
    // "not found".
    bool unlink(const std::string& job_id, std::string& error_msg);

    // Walk the entire root tree and remove every *.bin.zst.tmp file
    // (orphans from a crash or hard kill).  Returns the count of files
    // removed; sets error_msg on iteration failure.
    int sweep_tmp(std::string& error_msg);

    // Enumerate every job_id whose <ab>/<job_id>.bin.zst file currently
    // exists under <root>.  Used by JobHousekeeper's orphan sweep to
    // find files whose SQLite row was already deleted.
    std::vector<std::string> list_job_ids(std::string& error_msg) const;

    int level() const { return level_; }
    const std::string& root() const { return root_; }

private:
    std::string root_;
    int         level_;

    static std::string shard_for(const std::string& job_id);
};

} // namespace ikafssn
