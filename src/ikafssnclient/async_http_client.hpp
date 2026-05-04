#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "ikafssnclient/http_client.hpp"
#include "protocol/messages.hpp"

namespace ikafssn {

// Status snapshot returned by GET /api/v1/jobs/<id>.
struct AsyncJobStatus {
    std::string status;        // "queued" | "running" | "done" | "failed"
    std::string error_message;
    std::string fail_reason;
    int32_t     attempts = 0;
    int64_t     submitted_at = 0;
    int64_t     started_at = 0;
    int64_t     completed_at = 0;
};

// HTTP outcome categories for the new async REST endpoints.  Lets the
// poll loop tell apart "transient transport error" (retry next tick)
// from "remote says 404" (treat as terminally failed).
enum class AsyncHttpOutcome {
    kOk,
    kTransport,    // network failure, retry with backoff
    kNotFound,     // HTTP 404 (job missing — likely TTL expired)
    kConflict,     // HTTP 409 (duplicate job_id, or status != done)
    kBadRequest,   // HTTP 400 (server rejected payload)
    kServerError,  // HTTP 5xx
};

// POST /api/v1/jobs.  Caller supplies job_id (UUIDv4) — the body
// includes it so retries are idempotent: a successful resubmit returns
// 202 again, a duplicate returns 409 (treat as success — the job
// already exists).
AsyncHttpOutcome http_submit_job(const std::string& base_url,
                                 const std::string& job_id,
                                 const SearchRequest& req,
                                 std::string& error_msg,
                                 const HttpAuthConfig& auth = {});

// GET /api/v1/jobs/<id>.
AsyncHttpOutcome http_get_job_status(const std::string& base_url,
                                     const std::string& job_id,
                                     AsyncJobStatus& out,
                                     std::string& error_msg,
                                     const HttpAuthConfig& auth = {});

// GET /api/v1/jobs/<id>/result.  Returns the binary blob for
// deserialize(SearchResponse).
AsyncHttpOutcome http_get_job_result(const std::string& base_url,
                                     const std::string& job_id,
                                     std::vector<uint8_t>& out,
                                     std::string& error_msg,
                                     const HttpAuthConfig& auth = {});

} // namespace ikafssn
