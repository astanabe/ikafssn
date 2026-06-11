#pragma once

#include <atomic>
#include <string>
#include <vector>

#include "ikafssnclient/async_http_client.hpp"
#include "ikafssnclient/job_dir.hpp"
#include "io/result_writer.hpp"
#include "util/logger.hpp"

namespace ikafssn {

// Poll a group's outstanding jobs to completion, downloading their
// result blobs into `<job_id>.result.bin.zst` as each one finishes.
//
// Schedule: 30 seconds for the first 10 polls, then 60 seconds for
// subsequent polls (per design doc).  SIGINT (Ctrl+C) should toggle
// `g_sigint_flag` (set by the binary's signal handler) to make the
// loop return early; the persisted state in `<group_id>/` allows a
// later `-resume` to pick up where this run left off.
class PollLoop {
public:
    PollLoop(const std::string& jobs_root,
             const std::string& group_id,
             const HttpAuthConfig& auth,
             Logger& logger);

    // Block until every job in the group is `done` or `failed`, or
    // the SIGINT flag fires.  Returns true on full completion, false
    // when interrupted (caller exits 130) or on a fatal error.
    bool run();

    static std::atomic<bool>& sigint_flag();

private:
    std::string     root_;
    std::string     group_id_;
    HttpAuthConfig  auth_;
    Logger&         logger_;

    GroupMeta       group_;

    bool poll_once_();
    bool refresh_one_(ClientJobMeta& jm);
    int  next_delay_seconds_(int poll_count) const;
};

} // namespace ikafssn
