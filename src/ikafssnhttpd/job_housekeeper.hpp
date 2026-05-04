#pragma once

#include <atomic>
#include <condition_variable>
#include <mutex>
#include <thread>

#include "ikafssnhttpd/job_store.hpp"
#include "util/logger.hpp"

namespace ikafssn {

// Periodic TTL sweep over the job store.  Walks `delete_expired` at
// `interval = max(1, min(60, retention_time / 100))` seconds so that
// long retentions do not waste CPU and short retentions still tick at
// least once per second.
class JobHousekeeper {
public:
    JobHousekeeper(JobStore& store, Logger& logger);
    ~JobHousekeeper();

    JobHousekeeper(const JobHousekeeper&) = delete;
    JobHousekeeper& operator=(const JobHousekeeper&) = delete;

    void start(int retention_time_seconds);
    void stop();

private:
    JobStore& store_;
    Logger&   logger_;
    int       retention_time_seconds_ = 0;
    std::thread thread_;
    std::atomic<bool> stop_flag_{false};
    std::mutex mu_;
    std::condition_variable cv_;

    void loop_(int interval_seconds);
};

} // namespace ikafssn
