#include "ikafssnhttpd/job_housekeeper.hpp"

#include <algorithm>
#include <chrono>

namespace ikafssn {

JobHousekeeper::JobHousekeeper(JobStore& store,
                               QueryStore& queries,
                               ResultStore& results,
                               Logger& logger)
    : store_(store), queries_(queries), results_(results), logger_(logger) {}

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
    int tick = 0;
    while (!stop_flag_.load()) {
        {
            std::unique_lock<std::mutex> lk(mu_);
            cv_.wait_for(lk, std::chrono::seconds(interval_seconds),
                         [this] { return stop_flag_.load(); });
        }
        if (stop_flag_.load()) break;

        std::vector<std::string> deleted_ids;
        std::string err;
        int64_t n = store_.delete_expired(retention_time_seconds_,
                                          deleted_ids, err);
        if (!err.empty()) {
            logger_.error("Housekeeper delete_expired failed: %s", err.c_str());
        } else if (n > 0) {
            logger_.info("Housekeeper: removed %lld expired job(s)",
                         static_cast<long long>(n));
        }
        // Remove the matching result and (best-effort) query files.
        // The query file is normally already gone — JobWorker unlinks
        // it on mark_done / terminal mark_failed — but we still try
        // here as insurance against a worker crash between mark_done
        // and unlink.
        for (const auto& jid : deleted_ids) {
            std::string uerr;
            if (!results_.unlink(jid, uerr)) {
                logger_.error("Housekeeper: unlink result for %s: %s",
                              jid.c_str(), uerr.c_str());
            }
            std::string qerr;
            if (!queries_.unlink(jid, qerr)) {
                logger_.error("Housekeeper: unlink query for %s: %s",
                              jid.c_str(), qerr.c_str());
            }
        }

        if (++tick % kOrphanSweepEveryNTicks == 0) {
            orphan_sweep_();
        }
    }
}

void JobHousekeeper::orphan_sweep_() {
    // Result store sweep.
    {
        std::string err;
        auto ids = results_.list_job_ids(err);
        if (!err.empty()) {
            logger_.error("Housekeeper orphan sweep (results): list_job_ids: %s",
                          err.c_str());
        } else {
            int removed = 0;
            for (const auto& jid : ids) {
                JobMeta meta;
                std::string gerr;
                if (store_.get_status(jid, meta, gerr)) continue;
                if (!gerr.empty()) {
                    logger_.error("Housekeeper orphan sweep (results): "
                                  "get_status %s: %s",
                                  jid.c_str(), gerr.c_str());
                    continue;
                }
                std::string uerr;
                if (results_.unlink(jid, uerr)) {
                    removed++;
                } else {
                    logger_.error("Housekeeper orphan sweep (results): "
                                  "unlink %s: %s",
                                  jid.c_str(), uerr.c_str());
                }
            }
            if (removed > 0) {
                logger_.info("Housekeeper: orphan sweep removed %d "
                             "result file(s)", removed);
            }
        }
    }

    // Query store sweep.  An entry whose SQLite row reports 'done' or
    // 'failed' is also orphaned (the worker should have unlinked it on
    // terminal status, but a crash between mark_done and unlink leaves
    // the file behind).
    {
        std::string err;
        auto ids = queries_.list_job_ids(err);
        if (!err.empty()) {
            logger_.error("Housekeeper orphan sweep (queries): list_job_ids: %s",
                          err.c_str());
            return;
        }
        int removed = 0;
        for (const auto& jid : ids) {
            JobMeta meta;
            std::string gerr;
            bool has_row = store_.get_status(jid, meta, gerr);
            if (!gerr.empty()) {
                logger_.error("Housekeeper orphan sweep (queries): "
                              "get_status %s: %s",
                              jid.c_str(), gerr.c_str());
                continue;
            }
            bool should_remove = !has_row
                              || meta.status == JobStatus::kDone
                              || meta.status == JobStatus::kFailed;
            if (!should_remove) continue;
            std::string uerr;
            if (queries_.unlink(jid, uerr)) {
                removed++;
            } else {
                logger_.error("Housekeeper orphan sweep (queries): "
                              "unlink %s: %s",
                              jid.c_str(), uerr.c_str());
            }
        }
        if (removed > 0) {
            logger_.info("Housekeeper: orphan sweep removed %d "
                         "query file(s)", removed);
        }
    }
}

} // namespace ikafssn
