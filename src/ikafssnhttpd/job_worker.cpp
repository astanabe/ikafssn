#include "ikafssnhttpd/job_worker.hpp"

#include <chrono>
#include <utility>

#include "protocol/serializer.hpp"

namespace ikafssn {

JobWorker::JobWorker(JobStore& store,
                     ResultStore& results,
                     std::shared_ptr<BackendManager> manager,
                     Logger& logger,
                     int max_nretry)
    : store_(store)
    , results_(results)
    , manager_(std::move(manager))
    , logger_(logger)
    , max_nretry_(max_nretry > 0 ? max_nretry : 1) {}

JobWorker::~JobWorker() { stop(); }

void JobWorker::start(int n_workers) {
    if (!threads_.empty()) return;
    if (n_workers < 1) n_workers = 1;
    stop_flag_.store(false);
    threads_.reserve(static_cast<size_t>(n_workers));
    for (int i = 0; i < n_workers; i++) {
        threads_.emplace_back([this] { worker_loop_(); });
    }
}

void JobWorker::notify() { cv_.notify_all(); }

void JobWorker::stop() {
    if (threads_.empty()) return;
    stop_flag_.store(true);
    cv_.notify_all();
    for (auto& t : threads_) {
        if (t.joinable()) t.join();
    }
    threads_.clear();
}

void JobWorker::worker_loop_() {
    while (!stop_flag_.load()) {
        JobMeta meta;
        std::vector<uint8_t> request_blob;
        std::string err;
        bool got = store_.fetch_one_queued(meta, request_blob, err);
        if (!got) {
            if (!err.empty()) {
                logger_.error("JobWorker fetch failed: %s", err.c_str());
            }
            // Sleep until something pokes us, but cap at 5s so a missed
            // wake-up doesn't strand the queue forever.
            std::unique_lock<std::mutex> lk(cv_mu_);
            cv_.wait_for(lk, std::chrono::seconds(5),
                         [this] { return stop_flag_.load(); });
            continue;
        }

        process_one_(meta, request_blob);
    }
}

void JobWorker::process_one_(JobMeta& meta,
                             std::vector<uint8_t>& request_blob) {
    SearchRequest req;
    if (!deserialize(request_blob, req)) {
        std::string err;
        store_.mark_failed(meta.job_id,
                           "request_blob deserialize failed",
                           "request_blob_corrupt",
                           err);
        logger_.error("Job %s: request_blob deserialise failed",
                      meta.job_id.c_str());
        return;
    }

    SearchResponse resp;
    std::string backend_err;
    bool ok = manager_->route_search(req, resp, backend_err);

    if (ok) {
        auto blob = serialize(resp);
        std::string write_err;
        if (!results_.write(meta.job_id, blob.data(), blob.size(),
                            write_err)) {
            // ENOSPC / permission error / etc.: switch to mark_failed
            // so the client sees a clean failure rather than waiting
            // forever for a result that will never appear.
            std::string fail_reason = "result_write_failed: " + write_err;
            std::string mf_err;
            if (!store_.mark_failed(meta.job_id, write_err, fail_reason,
                                    mf_err)) {
                logger_.error("Job %s: mark_failed (after result_write) "
                              "failed: %s",
                              meta.job_id.c_str(), mf_err.c_str());
            } else {
                logger_.error("Job %s: result write failed (%s); marked failed",
                              meta.job_id.c_str(), write_err.c_str());
            }
            return;
        }
        std::string err;
        if (!store_.mark_done(meta.job_id, err)) {
            // The file is now an orphan; the housekeeper's periodic
            // sweep will clean it up.  Log loudly because this is rare.
            logger_.error("Job %s: mark_done failed (file written): %s",
                          meta.job_id.c_str(), err.c_str());
        }
        return;
    }

    // Failure path.
    int next_attempts = meta.attempts + 1;
    bool can_retry =
        (next_attempts < max_nretry_) && manager_->any_healthy();

    if (can_retry) {
        std::string err;
        if (!store_.requeue_for_retry(meta.job_id, backend_err, err)) {
            logger_.error("Job %s: requeue_for_retry failed: %s",
                          meta.job_id.c_str(), err.c_str());
        } else {
            logger_.info("Job %s: re-queued (attempts=%d, last_err=%s)",
                         meta.job_id.c_str(), next_attempts,
                         backend_err.c_str());
            cv_.notify_one();
        }
        return;
    }

    std::string fail_reason;
    if (!manager_->any_healthy()) {
        fail_reason = "backend_unreachable: all backends excluded";
    } else {
        fail_reason = "backend_error: " + backend_err;
    }

    std::string err;
    if (!store_.mark_failed(meta.job_id, backend_err, fail_reason, err)) {
        logger_.error("Job %s: mark_failed failed: %s",
                      meta.job_id.c_str(), err.c_str());
    } else {
        logger_.info("Job %s: failed (%s)",
                     meta.job_id.c_str(), fail_reason.c_str());
    }
}

} // namespace ikafssn
