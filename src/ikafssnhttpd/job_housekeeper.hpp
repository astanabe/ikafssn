#pragma once

#include <atomic>
#include <condition_variable>
#include <mutex>
#include <thread>

#include "ikafssnhttpd/job_store.hpp"
#include "ikafssnhttpd/result_store.hpp"
#include "util/logger.hpp"

namespace ikafssn {

// Periodic TTL sweep over the job store.  Walks `delete_expired` at
// `interval = max(1, min(60, retention_time / 100))` seconds so that
// long retentions do not waste CPU and short retentions still tick at
// least once per second.  After the SQLite rows are deleted, the
// matching ResultStore files are unlinked (best effort: log only).
//
// In addition, every `kOrphanSweepEveryNTicks` cycles the housekeeper
// scans the ResultStore for files whose SQLite row has already gone
// (e.g. mark_done failed after the file was written) and removes them.
class JobHousekeeper {
public:
    JobHousekeeper(JobStore& store, ResultStore& results, Logger& logger);
    ~JobHousekeeper();

    JobHousekeeper(const JobHousekeeper&) = delete;
    JobHousekeeper& operator=(const JobHousekeeper&) = delete;

    void start(int retention_time_seconds);
    void stop();

private:
    JobStore&    store_;
    ResultStore& results_;
    Logger&      logger_;
    int          retention_time_seconds_ = 0;
    std::thread thread_;
    std::atomic<bool> stop_flag_{false};
    std::mutex mu_;
    std::condition_variable cv_;

    static constexpr int kOrphanSweepEveryNTicks = 60;

    void loop_(int interval_seconds);
    void orphan_sweep_();
};

} // namespace ikafssn
