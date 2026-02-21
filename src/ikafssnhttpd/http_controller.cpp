#include "ikafssnhttpd/http_controller.hpp"

#include <thread>

#include <drogon/HttpAppFramework.h>
#include <json/json.h>

namespace ikafssn {

HttpController::HttpController(std::shared_ptr<BackendClient> backend)
    : backend_(std::move(backend)) {}

void HttpController::register_routes(const std::string& path_prefix) {
    std::string prefix = path_prefix;
    if (!prefix.empty() && prefix.back() == '/') {
        prefix.pop_back();
    }

    auto self = this;

    drogon::app().registerHandler(
        prefix + "/api/v1/search",
        [self](const drogon::HttpRequestPtr& req,
               std::function<void(const drogon::HttpResponsePtr&)>&& callback) {
            self->search(req, std::move(callback));
        },
        {drogon::Post});

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

void HttpController::search(
    const drogon::HttpRequestPtr& req,
    std::function<void(const drogon::HttpResponsePtr&)>&& callback) {

    // Parse JSON body
    auto json = req->getJsonObject();
    if (!json) {
        callback(make_error_response(drogon::k400BadRequest,
                                     "Invalid or missing JSON body"));
        return;
    }

    // Build SearchRequest from JSON
    SearchRequest sreq;
    const auto& j = *json;

    sreq.k = static_cast<uint8_t>(j.get("k", 0).asUInt());
    sreq.min_score = static_cast<uint16_t>(j.get("min_score", 0).asUInt());
    sreq.max_gap = static_cast<uint16_t>(j.get("max_gap", 0).asUInt());
    if (j.isMember("max_freq_frac") && j["max_freq_frac"].isDouble()) {
        double frac = j["max_freq_frac"].asDouble();
        if (frac > 0 && frac < 1.0) {
            sreq.max_freq_frac_x10000 = static_cast<uint16_t>(frac * 10000.0);
        }
    } else {
        sreq.max_freq = j.get("max_freq", 0).asUInt();
    }
    sreq.min_diag_hits = static_cast<uint8_t>(j.get("min_diag_hits", 0).asUInt());
    sreq.stage1_topn = static_cast<uint16_t>(j.get("stage1_topn", 0).asUInt());
    sreq.min_stage1_score = static_cast<uint16_t>(j.get("min_stage1_score", 0).asUInt());
    sreq.num_results = static_cast<uint16_t>(j.get("num_results", 0).asUInt());
    sreq.mode = static_cast<uint8_t>(j.get("mode", 0).asUInt());
    sreq.stage1_score_type = static_cast<uint8_t>(j.get("stage1_score", 0).asUInt());
    sreq.sort_score = static_cast<uint8_t>(j.get("sort_score", 0).asUInt());
    sreq.accept_qdegen = static_cast<uint8_t>(j.get("accept_qdegen", 0).asUInt());
    sreq.strand = static_cast<int8_t>(j.get("strand", 0).asInt());

    // Seqidlist mode
    std::string mode_str = j.get("seqidlist_mode", "none").asString();
    if (mode_str == "include") {
        sreq.seqidlist_mode = SeqidlistMode::kInclude;
    } else if (mode_str == "exclude") {
        sreq.seqidlist_mode = SeqidlistMode::kExclude;
    } else {
        sreq.seqidlist_mode = SeqidlistMode::kNone;
    }

    // Seqids array
    if (j.isMember("seqids") && j["seqids"].isArray()) {
        for (const auto& s : j["seqids"]) {
            sreq.seqids.push_back(s.asString());
        }
    }

    // Queries array (required)
    if (!j.isMember("queries") || !j["queries"].isArray() ||
        j["queries"].empty()) {
        callback(make_error_response(drogon::k400BadRequest,
                                     "Missing or empty 'queries' array"));
        return;
    }

    for (const auto& q : j["queries"]) {
        if (!q.isMember("query_id") || !q.isMember("sequence")) {
            callback(make_error_response(
                drogon::k400BadRequest,
                "Each query must have 'query_id' and 'sequence' fields"));
            return;
        }
        QueryEntry entry;
        entry.query_id = q["query_id"].asString();
        entry.sequence = q["sequence"].asString();
        if (entry.sequence.empty()) {
            callback(make_error_response(drogon::k400BadRequest,
                                         "Query sequence must not be empty"));
            return;
        }
        sreq.queries.push_back(std::move(entry));
    }

    // Offload blocking backend I/O to a worker thread.
    // Drogon event loop threads must not be blocked.
    auto backend = backend_;
    auto cb = std::make_shared<std::function<void(const drogon::HttpResponsePtr&)>>(
        std::move(callback));

    std::thread([backend, sreq = std::move(sreq), cb]() {
        SearchResponse sresp;
        std::string error_msg;
        if (!backend->search(sreq, sresp, error_msg)) {
            (*cb)(make_error_response(drogon::k502BadGateway, error_msg));
            return;
        }

        // Build JSON response
        Json::Value result;
        result["status"] = (sresp.status == 0) ? "success" : "error";
        result["k"] = sresp.k;
        result["mode"] = sresp.mode;
        result["stage1_score_type"] = sresp.stage1_score_type;
        const char* s1name = (sresp.stage1_score_type == 2) ? "matchscore" : "coverscore";

        Json::Value results_arr(Json::arrayValue);
        for (const auto& qr : sresp.results) {
            Json::Value qobj;
            qobj["query_id"] = qr.query_id;

            Json::Value hits_arr(Json::arrayValue);
            for (const auto& hit : qr.hits) {
                Json::Value hobj;
                hobj["accession"] = hit.accession;
                hobj["strand"] = (hit.strand == 0) ? "+" : "-";
                if (sresp.mode != 1) {
                    hobj["q_start"] = hit.q_start;
                    hobj["q_end"] = hit.q_end;
                    hobj["s_start"] = hit.s_start;
                    hobj["s_end"] = hit.s_end;
                }
                hobj[s1name] = hit.stage1_score;
                if (sresp.mode != 1) {
                    hobj["chainscore"] = hit.score;
                }
                hobj["volume"] = hit.volume;
                hits_arr.append(std::move(hobj));
            }
            qobj["hits"] = std::move(hits_arr);
            if (qr.skipped != 0) {
                qobj["skipped"] = true;
            }
            results_arr.append(std::move(qobj));
        }
        result["results"] = std::move(results_arr);

        if (!sresp.rejected_query_ids.empty()) {
            Json::Value rejected_arr(Json::arrayValue);
            for (const auto& qid : sresp.rejected_query_ids)
                rejected_arr.append(qid);
            result["rejected_query_ids"] = std::move(rejected_arr);
        }

        auto resp = drogon::HttpResponse::newHttpJsonResponse(std::move(result));
        (*cb)(resp);
    }).detach();
}

void HttpController::health(
    const drogon::HttpRequestPtr& req,
    std::function<void(const drogon::HttpResponsePtr&)>&& callback) {

    auto backend = backend_;
    auto cb = std::make_shared<std::function<void(const drogon::HttpResponsePtr&)>>(
        std::move(callback));

    std::thread([backend, cb]() {
        HealthResponse hresp;
        std::string error_msg;
        if (!backend->health_check(hresp, error_msg)) {
            (*cb)(make_error_response(drogon::k502BadGateway, error_msg));
            return;
        }

        Json::Value result;
        result["status"] = (hresp.status == 0) ? "ok" : "error";

        auto resp = drogon::HttpResponse::newHttpJsonResponse(std::move(result));
        (*cb)(resp);
    }).detach();
}

void HttpController::info(
    const drogon::HttpRequestPtr& req,
    std::function<void(const drogon::HttpResponsePtr&)>&& callback) {

    auto backend = backend_;
    auto cb = std::make_shared<std::function<void(const drogon::HttpResponsePtr&)>>(
        std::move(callback));

    std::thread([backend, cb]() {
        InfoResponse iresp;
        std::string error_msg;
        if (!backend->info(iresp, error_msg)) {
            (*cb)(make_error_response(drogon::k502BadGateway, error_msg));
            return;
        }

        Json::Value result;
        result["status"] = (iresp.status == 0) ? "success" : "error";
        result["default_k"] = iresp.default_k;

        Json::Value groups_arr(Json::arrayValue);
        for (const auto& g : iresp.groups) {
            Json::Value gobj;
            gobj["k"] = g.k;
            gobj["kmer_type"] = (g.kmer_type == 0) ? "uint16" : "uint32";

            uint64_t group_total_sequences = 0;
            uint64_t group_total_postings = 0;

            Json::Value vols_arr(Json::arrayValue);
            for (const auto& v : g.volumes) {
                Json::Value vobj;
                vobj["volume_index"] = v.volume_index;
                vobj["num_sequences"] = v.num_sequences;
                vobj["total_postings"] = static_cast<Json::UInt64>(v.total_postings);
                vobj["db_name"] = v.db_name;
                vols_arr.append(std::move(vobj));

                group_total_sequences += v.num_sequences;
                group_total_postings += v.total_postings;
            }
            gobj["volumes"] = std::move(vols_arr);
            gobj["num_volumes"] = static_cast<Json::UInt>(g.volumes.size());
            gobj["total_sequences"] = static_cast<Json::UInt64>(group_total_sequences);
            gobj["total_postings"] = static_cast<Json::UInt64>(group_total_postings);
            groups_arr.append(std::move(gobj));
        }
        result["kmer_groups"] = std::move(groups_arr);

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
