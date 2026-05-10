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
#include "ikafssnhttpd/query_store.hpp"
#include "ikafssnhttpd/result_store.hpp"
#include "protocol/messages.hpp"

namespace ikafssn {

// HTTP REST API controller.  Job-broker model:
//
//   POST   /api/v1/jobs              submit (returns 202 + {"job_id"})
//   GET    /api/v1/jobs/<id>          poll status JSON
//   GET    /api/v1/jobs/<id>/result   download serialized SearchResponse
//
// The controller owns no worker pool of its own; every POST validates
// the request and persists the body to QueryStore, then inserts a row
// into JobStore and notifies the (separately-owned) JobWorker pool.
class HttpController {
public:
    HttpController(std::shared_ptr<BackendManager> manager,
                   JobStore& store,
                   JobWorker& worker,
                   QueryStore& queries,
                   ResultStore& results);

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

private:
    std::shared_ptr<BackendManager> manager_;
    JobStore*    store_;
    JobWorker*   worker_;
    QueryStore*  queries_;
    ResultStore* results_;

    static drogon::HttpResponsePtr make_error_response(
        drogon::HttpStatusCode status, const std::string& message);
};

} // namespace ikafssn
