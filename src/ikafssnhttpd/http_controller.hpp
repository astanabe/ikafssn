#pragma once

#include <functional>
#include <memory>
#include <string>

#include <drogon/HttpRequest.h>
#include <drogon/HttpResponse.h>
#include <drogon/HttpTypes.h>

#include "ikafssnhttpd/backend_client.hpp"

namespace ikafssn {

// HTTP REST API controller.
// Translates HTTP JSON requests to binary protocol messages (via BackendClient)
// and binary protocol responses back to HTTP JSON responses.
class HttpController {
public:
    explicit HttpController(std::shared_ptr<BackendClient> backend);

    // Register HTTP routes with Drogon. Must be called before app().run().
    void register_routes(const std::string& path_prefix);

    // POST /api/v1/search
    void search(const drogon::HttpRequestPtr& req,
                std::function<void(const drogon::HttpResponsePtr&)>&& callback);

    // GET /api/v1/health
    void health(const drogon::HttpRequestPtr& req,
                std::function<void(const drogon::HttpResponsePtr&)>&& callback);

    // GET /api/v1/info
    void info(const drogon::HttpRequestPtr& req,
              std::function<void(const drogon::HttpResponsePtr&)>&& callback);

private:
    std::shared_ptr<BackendClient> backend_;

    static drogon::HttpResponsePtr make_error_response(
        drogon::HttpStatusCode status, const std::string& message);
};

} // namespace ikafssn
