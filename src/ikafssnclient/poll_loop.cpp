#include "ikafssnclient/poll_loop.hpp"

#include <chrono>
#include <thread>

namespace ikafssn {

namespace {

int64_t now_unix() {
    using namespace std::chrono;
    return duration_cast<seconds>(system_clock::now().time_since_epoch()).count();
}

bool is_terminal(const std::string& status) {
    return status == "done" || status == "failed" || status == "timeout";
}

} // namespace

std::atomic<bool>& PollLoop::sigint_flag() {
    static std::atomic<bool> flag{false};
    return flag;
}

PollLoop::PollLoop(const std::string& jobs_root,
                   const std::string& group_id,
                   const HttpAuthConfig& auth,
                   Logger& logger)
    : root_(jobs_root)
    , group_id_(group_id)
    , auth_(auth)
    , logger_(logger) {}

int PollLoop::next_delay_seconds_(int poll_count) const {
    return (poll_count < 10) ? 30 : 60;
}

bool PollLoop::run() {
    std::string err;
    if (!read_group_meta(root_, group_id_, group_, err)) {
        logger_.error("PollLoop: cannot read group meta: %s", err.c_str());
        return false;
    }

    int poll_count = 0;
    while (!sigint_flag().load()) {
        bool all_done = poll_once_();
        if (all_done) return true;

        int delay = next_delay_seconds_(poll_count++);
        for (int i = 0; i < delay && !sigint_flag().load(); i++) {
            std::this_thread::sleep_for(std::chrono::seconds(1));
        }
    }
    logger_.info("PollLoop: SIGINT, leaving %s in place for -resume",
                 (root_ + "/" + group_id_).c_str());
    return false;
}

bool PollLoop::poll_once_() {
    int n_terminal = 0;
    int n_total    = 0;
    int cnt_q = 0, cnt_r = 0, cnt_d = 0, cnt_f = 0;

    for (const auto& job_id : group_.job_ids) {
        n_total++;
        JobMeta jm;
        std::string err;
        if (!read_job_meta(root_, group_id_, job_id, jm, err)) {
            logger_.error("PollLoop: read_job_meta(%s): %s",
                          job_id.c_str(), err.c_str());
            continue;
        }
        if (!is_terminal(jm.status)) {
            refresh_one_(jm);
        }

        if (jm.status == "queued")       cnt_q++;
        else if (jm.status == "running") cnt_r++;
        else if (jm.status == "done")    cnt_d++;
        else if (jm.status == "failed" || jm.status == "timeout") cnt_f++;

        if (is_terminal(jm.status)) n_terminal++;
    }

    group_.cnt_queued  = cnt_q;
    group_.cnt_running = cnt_r;
    group_.cnt_done    = cnt_d;
    group_.cnt_failed  = cnt_f;
    {
        std::string err;
        if (!write_group_meta(root_, group_, err)) {
            logger_.error("PollLoop: update group_meta: %s", err.c_str());
        }
    }

    logger_.info("Poll: queued=%d running=%d done=%d failed=%d (of %d)",
                 cnt_q, cnt_r, cnt_d, cnt_f, n_total);
    return n_terminal == n_total;
}

bool PollLoop::refresh_one_(JobMeta& jm) {
    AsyncJobStatus st;
    std::string err;
    auto outcome = http_get_job_status(group_.httpd_url, jm.job_id,
                                        st, err, auth_);
    if (outcome == AsyncHttpOutcome::kNotFound) {
        // Server forgot it — treat as TTL timeout.
        jm.status = "timeout";
        jm.fail_reason = "timeout: exceeded retention_time";
        jm.completed_at = now_unix();
        std::string werr;
        write_job_meta(root_, jm, werr);
        return true;
    }
    if (outcome != AsyncHttpOutcome::kOk) {
        jm.last_polled_at = now_unix();
        std::string werr;
        write_job_meta(root_, jm, werr);
        logger_.debug("Poll %s: transport %s", jm.job_id.c_str(), err.c_str());
        return false;
    }

    jm.status        = st.status;
    jm.attempts      = st.attempts;
    jm.error_message = st.error_message;
    jm.fail_reason   = st.fail_reason;
    jm.last_polled_at = now_unix();
    if (st.completed_at) jm.completed_at = st.completed_at;

    if (st.status == "done") {
        std::vector<uint8_t> blob;
        std::string rerr;
        auto rc = http_get_job_result(group_.httpd_url, jm.job_id,
                                       blob, rerr, auth_);
        if (rc == AsyncHttpOutcome::kOk) {
            std::string werr;
            if (!write_job_result(root_, group_id_, jm.job_id, blob, werr)) {
                logger_.error("Poll: cache result for %s: %s",
                              jm.job_id.c_str(), werr.c_str());
            }
        } else {
            logger_.error("Poll: GET result %s failed: %s",
                          jm.job_id.c_str(), rerr.c_str());
        }
    }
    std::string werr;
    write_job_meta(root_, jm, werr);
    return is_terminal(jm.status);
}

} // namespace ikafssn
