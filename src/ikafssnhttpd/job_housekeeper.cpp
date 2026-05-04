#include "ikafssnhttpd/job_housekeeper.hpp"

#include <algorithm>
#include <chrono>

namespace ikafssn {

JobHousekeeper::JobHousekeeper(JobStore& store, Logger& logger)
    : store_(store), logger_(logger) {}

JobHousekeeper::~JobHousekeeper() { stop(); }

void JobHousekeeper::start(int retention_time_seconds) {
    if (thread_.joinable()) return;
    retention_time_seconds_ = retention_time_seconds;
    int interval = std::max(1,
                            std::min(60, retention_time_seconds_ / 100));
    if (interval <= 0) interval = 1;
    stop_flag_.store(false);
    thread_ = std::thread([this, interval] { loop_(interval); });
}

void JobHousekeeper::stop() {
    stop_flag_.store(true);
    cv_.notify_all();
    if (thread_.joinable()) thread_.join();
}

void JobHousekeeper::loop_(int interval_seconds) {
    while (!stop_flag_.load()) {
        {
            std::unique_lock<std::mutex> lk(mu_);
            cv_.wait_for(lk, std::chrono::seconds(interval_seconds),
                         [this] { return stop_flag_.load(); });
        }
        if (stop_flag_.load()) break;

        std::string err;
        int64_t n = store_.delete_expired(retention_time_seconds_, err);
        if (!err.empty()) {
            logger_.error("Housekeeper delete_expired failed: %s", err.c_str());
        } else if (n > 0) {
            logger_.info("Housekeeper: removed %lld expired job(s)",
                         static_cast<long long>(n));
        }
    }
}

} // namespace ikafssn
