#pragma once

#include <string>

#include "protocol/messages.hpp"

namespace ikafssn {

// Parse the JSON body of POST /api/v1/jobs into (job_id, SearchRequest).
//
// The same parser is used by HttpController::submit_job (validating the
// upload before queuing) and by JobWorker::process_one_ (re-parsing the
// stored query before dispatch).  Lives in the protocol library so it
// can be linked from both the controller and the worker.
//
// On failure, `error_msg` is set to a human-readable reason and the
// function returns false.  The contents of `job_id` / `req` are
// undefined on failure.
bool parse_search_request_json(const std::string& body,
                               std::string& job_id,
                               SearchRequest& req,
                               std::string& error_msg);

} // namespace ikafssn
