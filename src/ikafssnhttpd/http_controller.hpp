#pragma once

#include <functional>
#include <memory>
#include <string>

#include <drogon/HttpRequest.h>
#include <drogon/HttpResponse.h>
#include <drogon/HttpTypes.h>

#include "ikafssnhttpd/backend_manager.hpp"
#include "ikafssnhttpd/job_store.hpp"
#include "ikafssnhttpd/job_worker.hpp"
#include "protocol/messages.hpp"

namespace ikafssn {

// HTTP REST API controller.  v9 of the API replaces the synchronous
// POST /api/v1/search with a job-broker model:
//
//   POST   /api/v1/jobs              submit (returns 202 + {"job_id"})
//   GET    /api/v1/jobs/<id>          poll status JSON
//   GET    /api/v1/jobs/<id>/result   download serialized SearchResponse
//
// The controller owns no worker pool of its own; every POST simply
// inserts a row into JobStore and notifies the (separately-owned)
// JobWorker pool.
class HttpController {
public:
    HttpController(std::shared_ptr<BackendManager> manager,
                   JobStore& store,
                   JobWorker& worker);

    // Register HTTP routes with Drogon.  Must be called before app().run().
    void register_routes(const std::string& path_prefix);

    void submit_job(const drogon::HttpRequestPtr& req,
                    std::function<void(const drogon::HttpResponsePtr&)>&& callback);

    void get_job_status(const drogon::HttpRequestPtr& req,
                        std::function<void(const drogon::HttpResponsePtr&)>&& callback,
                        const std::string& job_id);

    void get_job_result(const drogon::HttpRequestPtr& req,
                        std::function<void(const drogon::HttpResponsePtr&)>&& callback,
                        const std::string& job_id);

    void health(const drogon::HttpRequestPtr& req,
                std::function<void(const drogon::HttpResponsePtr&)>&& callback);

    void info(const drogon::HttpRequestPtr& req,
              std::function<void(const drogon::HttpResponsePtr&)>&& callback);

    // Parse the JSON body of POST /api/v1/jobs into (job_id, SearchRequest).
    // On failure, sets `error_msg` and returns false.  Exposed for tests.
    static bool parse_search_request_json(const std::string& body,
                                          std::string& job_id,
                                          SearchRequest& req,
                                          std::string& error_msg);

private:
    std::shared_ptr<BackendManager> manager_;
    JobStore*  store_;
    JobWorker* worker_;

    static drogon::HttpResponsePtr make_error_response(
        drogon::HttpStatusCode status, const std::string& message);
};

} // namespace ikafssn
