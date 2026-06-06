#include "ikafssnhttpd/http_controller.hpp"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <climits>
#include <thread>

#include <drogon/HttpAppFramework.h>
#include <json/json.h>

#include "ikafssnhttpd/search_request_json.hpp"
#include "protocol/info_format.hpp"
#include "protocol/messages.hpp"
#include "util/zstd_oneshot.hpp"

namespace ikafssn {

HttpController::HttpController(std::shared_ptr<BackendManager> manager,
                               JobStore& store,
                               JobWorker& worker,
                               QueryStore& queries,
                               ResultStore& results)
    : manager_(std::move(manager))
    , store_(&store)
    , worker_(&worker)
    , queries_(&queries)
    , results_(&results) {}

namespace {

// Lower-cased media-type extracted from a Content-Type header (the part
// before the first ';').  Leading and trailing whitespace is trimmed.
std::string extract_media_type(const std::string& header) {
    size_t end = header.find(';');
    std::string mt = header.substr(0, end);
    auto not_space = [](unsigned char c) { return !std::isspace(c); };
    auto first = std::find_if(mt.begin(), mt.end(), not_space);
    auto last  = std::find_if(mt.rbegin(), mt.rend(), not_space).base();
    if (first >= last) return {};
    std::string out(first, last);
    for (auto& c : out) c = static_cast<char>(std::tolower(c));
    return out;
}

} // namespace

void HttpController::register_routes(const std::string& path_prefix) {
    std::string prefix = path_prefix;
    if (!prefix.empty() && prefix.back() == '/') {
        prefix.pop_back();
    }

    auto self = this;

    drogon::app().registerHandler(
        prefix + "/api/v1/jobs",
        [self](const drogon::HttpRequestPtr& req,
               std::function<void(const drogon::HttpResponsePtr&)>&& callback) {
            self->submit_job(req, std::move(callback));
        },
        {drogon::Post});

    drogon::app().registerHandler(
        prefix + "/api/v1/jobs/{job_id}",
        [self](const drogon::HttpRequestPtr& req,
               std::function<void(const drogon::HttpResponsePtr&)>&& callback,
               const std::string& job_id) {
            self->get_job_status(req, std::move(callback), job_id);
        },
        {drogon::Get});

    drogon::app().registerHandler(
        prefix + "/api/v1/jobs/{job_id}/result",
        [self](const drogon::HttpRequestPtr& req,
               std::function<void(const drogon::HttpResponsePtr&)>&& callback,
               const std::string& job_id) {
            self->get_job_result(req, std::move(callback), job_id);
        },
        {drogon::Get});

    drogon::app().registerHandler(
        prefix + "/api/v1/health",
        [self](const drogon::HttpRequestPtr& req,
               std::function<void(const drogon::HttpResponsePtr&)>&& callback) {
            self->health(req, std::move(callback));
        },
        {drogon::Get});

    drogon::app().registerHandler(
        prefix + "/api/v1/info",
        [self](const drogon::HttpRequestPtr& req,
               std::function<void(const drogon::HttpResponsePtr&)>&& callback) {
            self->info(req, std::move(callback));
        },
        {drogon::Get});
}

void HttpController::submit_job(
    const drogon::HttpRequestPtr& req,
    std::function<void(const drogon::HttpResponsePtr&)>&& callback) {

    std::string raw_body(req->getBody());
    std::string content_type = req->getHeader("Content-Type");
    std::string media_type = extract_media_type(content_type);

    // ikafssnclient always sends application/zstd.  application/json (and
    // an empty Content-Type) is also accepted so external tools — curl
    // and any one-off scripts — can submit a job without having to pipe
    // the body through `zstd` first.  This is NOT a back-compat path for
    // the previous ikafssnclient release; it is a documented entry point
    // for ad-hoc HTTP clients.
    std::string body;
    bool client_sent_zstd = false;
    if (media_type == "application/zstd") {
        std::vector<uint8_t> decompressed;
        std::string zerr;
        if (!zstd_decompress(raw_body.data(), raw_body.size(),
                             decompressed, zerr)) {
            callback(make_error_response(drogon::k400BadRequest,
                "invalid zstd frame: " + zerr));
            return;
        }
        body.assign(reinterpret_cast<const char*>(decompressed.data()),
                    decompressed.size());
        client_sent_zstd = true;
    } else if (media_type.empty() || media_type == "application/json") {
        body = raw_body;
    } else {
        callback(make_error_response(drogon::k415UnsupportedMediaType,
            "unsupported Content-Type: " + content_type
            + "; expected application/zstd or application/json"));
        return;
    }

    std::string job_id;
    SearchRequest sreq;
    std::string err;
    if (!parse_search_request_json(body, job_id, sreq, err)) {
        callback(make_error_response(drogon::k400BadRequest, err));
        return;
    }

    // Validate db / k / mode against merged backend capabilities.
    {
        auto merged = manager_->merged_info();
        std::string verr = validate_info(merged, sreq.db,
                                          sreq.k, sreq.mode, false,
                                          sreq.t, sreq.template_type);
        if (!verr.empty()) {
            callback(make_error_response(drogon::k400BadRequest, verr));
            return;
        }
    }

    // Persist the request body to QueryStore *before* the SQLite INSERT
    // so a worker that wakes up on the new row always finds the file.
    // If INSERT fails (duplicate, sqlite error), unlink the file we
    // just wrote so the queue does not accumulate orphans.
    {
        std::string write_err;
        bool wrote;
        if (client_sent_zstd) {
            wrote = queries_->write_compressed_passthrough(
                job_id, raw_body.data(), raw_body.size(), write_err);
        } else {
            wrote = queries_->write_plain(
                job_id, body.data(), body.size(), write_err);
        }
        if (!wrote) {
            callback(make_error_response(drogon::k500InternalServerError,
                "Failed to persist query body: " + write_err));
            return;
        }
    }

    int64_t submitted_at = 0;
    bool duplicate = false;
    std::string store_err;
    if (!store_->insert_job(job_id, sreq.db,
                            static_cast<int32_t>(sreq.queries.size()),
                            submitted_at, duplicate, store_err)) {
        std::string ignored;
        queries_->unlink(job_id, ignored);
        if (duplicate) {
            callback(make_error_response(drogon::k409Conflict,
                "Job already exists: " + job_id));
        } else {
            callback(make_error_response(drogon::k500InternalServerError,
                "Failed to insert job: " + store_err));
        }
        return;
    }
    worker_->notify();

    Json::Value ok;
    ok["job_id"] = job_id;
    ok["status"] = "queued";
    ok["submitted_at"] = static_cast<Json::Int64>(submitted_at);
    auto resp = drogon::HttpResponse::newHttpJsonResponse(std::move(ok));
    resp->setStatusCode(drogon::k202Accepted);
    callback(resp);
}

void HttpController::get_job_status(
    const drogon::HttpRequestPtr&,
    std::function<void(const drogon::HttpResponsePtr&)>&& callback,
    const std::string& job_id) {

    JobMeta meta;
    std::string err;
    if (!store_->get_status(job_id, meta, err)) {
        if (err.empty()) {
            callback(make_error_response(drogon::k404NotFound,
                "Unknown job_id: " + job_id));
        } else {
            callback(make_error_response(drogon::k500InternalServerError, err));
        }
        return;
    }
    Json::Value out;
    out["job_id"] = meta.job_id;
    out["status"] = job_status_str(meta.status);
    out["attempts"] = meta.attempts;
    out["submitted_at"] = static_cast<Json::Int64>(meta.submitted_at);
    if (meta.started_at)   out["started_at"]   = static_cast<Json::Int64>(meta.started_at);
    if (meta.completed_at) out["completed_at"] = static_cast<Json::Int64>(meta.completed_at);
    if (!meta.error_message.empty()) out["error_message"] = meta.error_message;
    if (!meta.fail_reason.empty())   out["fail_reason"]   = meta.fail_reason;
    if (!meta.db.empty())            out["db"]            = meta.db;
    out["n_seqs"] = meta.n_seqs;

    auto resp = drogon::HttpResponse::newHttpJsonResponse(std::move(out));
    callback(resp);
}

void HttpController::get_job_result(
    const drogon::HttpRequestPtr&,
    std::function<void(const drogon::HttpResponsePtr&)>&& callback,
    const std::string& job_id) {

    JobMeta meta;
    std::string err;
    if (!store_->get_status(job_id, meta, err)) {
        if (err.empty()) {
            callback(make_error_response(drogon::k404NotFound,
                "Unknown job_id: " + job_id));
        } else {
            callback(make_error_response(drogon::k500InternalServerError, err));
        }
        return;
    }
    if (meta.status != JobStatus::kDone) {
        callback(make_error_response(drogon::k409Conflict,
            "Job not yet complete"));
        return;
    }
    if (!results_->exists(job_id)) {
        // Race with the housekeeper: the SQLite row still says 'done'
        // but the file has already been swept (or never landed because
        // mark_done committed before the file was fsynced — extremely
        // unlikely given the worker order, but possible if an operator
        // hand-deleted the file).
        callback(make_error_response(drogon::k404NotFound,
            "result file missing for job " + job_id
            + ": housekeeper raced"));
        return;
    }

    auto resp = drogon::HttpResponse::newFileResponse(
        results_->path_for(job_id),
        /*attachmentFileName=*/"",
        drogon::CT_CUSTOM,
        /*typeString=*/"application/zstd");
    callback(resp);
}

void HttpController::health(
    const drogon::HttpRequestPtr&,
    std::function<void(const drogon::HttpResponsePtr&)>&& callback) {

    auto manager = manager_;
    auto cb = std::make_shared<std::function<void(const drogon::HttpResponsePtr&)>>(
        std::move(callback));
    std::thread([manager, cb]() {
        std::string error_msg;
        if (!manager->check_any_health(error_msg)) {
            (*cb)(make_error_response(drogon::k502BadGateway, error_msg));
            return;
        }
        Json::Value result;
        result["status"] = "ok";
        auto resp = drogon::HttpResponse::newHttpJsonResponse(std::move(result));
        (*cb)(resp);
    }).detach();
}

void HttpController::info(
    const drogon::HttpRequestPtr&,
    std::function<void(const drogon::HttpResponsePtr&)>&& callback) {

    auto manager = manager_;
    auto cb = std::make_shared<std::function<void(const drogon::HttpResponsePtr&)>>(
        std::move(callback));
    std::thread([manager, cb]() {
        auto result = manager->build_info_json();
        auto resp = drogon::HttpResponse::newHttpJsonResponse(std::move(result));
        (*cb)(resp);
    }).detach();
}

drogon::HttpResponsePtr HttpController::make_error_response(
    drogon::HttpStatusCode status, const std::string& message) {
    Json::Value body;
    body["error"] = message;
    auto resp = drogon::HttpResponse::newHttpJsonResponse(std::move(body));
    resp->setStatusCode(status);
    return resp;
}

} // namespace ikafssn
