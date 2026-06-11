#pragma once

#include <atomic>
#include <condition_variable>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "ikafssnhttpd/backend_manager.hpp"
#include "ikafssnhttpd/job_store.hpp"
#include "ikafssnhttpd/query_store.hpp"
#include "ikafssnhttpd/result_store.hpp"
#include "util/logger.hpp"

namespace ikafssn {

// Background worker pool that drains the JobStore queue.  Each worker
// fetches one queued job, reads its JSON body from the QueryStore,
// re-parses it into a SearchRequest, hands the request to
// BackendManager::route_search, and either marks the job as done,
// requeues it for retry, or marks it terminally failed.  After a
// terminal status (done or non-retry failed), the per-job query file
// is unlinked so the QueryStore only ever contains queued / running
// jobs.
//
// The retry policy lives entirely in this layer: BackendManager::route_search
// performs a single attempt, and retry counting is owned by attempts++
// in the JobStore.
class JobWorker {
public:
    JobWorker(JobStore& store,
              QueryStore& queries,
              ResultStore& results,
              std::shared_ptr<BackendManager> manager,
              Logger& logger,
              int max_nretry);

    ~JobWorker();

    JobWorker(const JobWorker&) = delete;
    JobWorker& operator=(const JobWorker&) = delete;

    // Spawn `n_workers` threads.  Idempotent: a second call with the
    // pool already running is a no-op.
    void start(int n_workers);

    // Wake any sleeping workers (e.g. after a new POST inserts a row).
    void notify();

    // Stop and join all workers.  Idempotent.
    void stop();

private:
    JobStore&                       store_;
    QueryStore&                     queries_;
    ResultStore&                    results_;
    std::shared_ptr<BackendManager> manager_;
    Logger&                         logger_;
    int                             max_nretry_;

    std::vector<std::thread>  threads_;
    std::atomic<bool>         stop_flag_{false};
    std::mutex                cv_mu_;
    std::condition_variable   cv_;

    void worker_loop_();
    void process_one_(HttpdJobMeta& meta);
};

} // namespace ikafssn
