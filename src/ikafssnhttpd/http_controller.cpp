#include "ikafssnhttpd/http_controller.hpp"

#include <chrono>
#include <climits>
#include <thread>

#include <drogon/HttpAppFramework.h>
#include <json/json.h>

#include "protocol/info_format.hpp"

namespace ikafssn {

HttpController::HttpController(std::shared_ptr<BackendClient> backend)
    : backend_(std::move(backend)) {}

bool HttpController::init_cache(int timeout_seconds) {
    auto start = std::chrono::steady_clock::now();
    int delay = 1;  // exponential backoff: 1, 2, 4, 8, ...

    while (true) {
        InfoResponse info;
        std::string error_msg;
        if (backend_->info(info, error_msg)) {
            std::lock_guard<std::mutex> lock(cache_mutex_);
            cached_info_ = std::move(info);
            cache_valid_ = true;
            return true;
        }

        auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(
            std::chrono::steady_clock::now() - start).count();
        if (elapsed >= timeout_seconds) {
            return false;
        }

        std::this_thread::sleep_for(std::chrono::seconds(delay));
        delay = std::min(delay * 2, 8);
    }
}

void HttpController::update_cache(const InfoResponse& info) {
    std::lock_guard<std::mutex> lock(cache_mutex_);
    cached_info_ = info;
    cache_valid_ = true;
}

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
    sreq.stage2_min_score = static_cast<uint16_t>(j.get("stage2_min_score", 0).asUInt());
    if (j.isMember("has_stage2_min_score") && j["has_stage2_min_score"].asBool()) {
        sreq.has_stage2_min_score = 1;
    }
    sreq.stage2_max_gap = static_cast<uint16_t>(j.get("stage2_max_gap", 0).asUInt());
    sreq.stage2_max_lookback = static_cast<uint16_t>(j.get("stage2_max_lookback", 0).asUInt());
    if (j.isMember("stage1_max_freq_frac") && j["stage1_max_freq_frac"].isDouble()) {
        double frac = j["stage1_max_freq_frac"].asDouble();
        if (frac > 0 && frac < 1.0) {
            sreq.stage1_max_freq_frac_x10000 = static_cast<uint16_t>(frac * 10000.0);
        }
    } else {
        sreq.stage1_max_freq = j.get("stage1_max_freq", 0).asUInt();
    }
    if (j.isMember("stage1_min_score_frac") && j["stage1_min_score_frac"].isDouble()) {
        double frac = j["stage1_min_score_frac"].asDouble();
        if (frac > 0 && frac < 1.0) {
            sreq.stage1_min_score_frac_x10000 = static_cast<uint16_t>(frac * 10000.0);
        }
    }
    sreq.stage2_min_diag_hits = static_cast<uint8_t>(j.get("stage2_min_diag_hits", 0).asUInt());
    sreq.stage1_topn = static_cast<uint16_t>(j.get("stage1_topn", 0).asUInt());
    sreq.stage1_min_score = static_cast<uint16_t>(j.get("stage1_min_score", 0).asUInt());
    sreq.num_results = static_cast<uint16_t>(j.get("num_results", 0).asUInt());
    sreq.mode = static_cast<uint8_t>(j.get("mode", 0).asUInt());
    sreq.stage1_score = static_cast<uint8_t>(j.get("stage1_score", 0).asUInt());
    sreq.accept_qdegen = static_cast<uint8_t>(j.get("accept_qdegen", 1).asUInt());
    sreq.strand = static_cast<int8_t>(j.get("strand", 0).asInt());

    // Stage 3 parameters
    sreq.stage3_traceback = static_cast<uint8_t>(j.get("stage3_traceback", 0).asUInt());
    sreq.stage3_gapopen = j.isMember("stage3_gapopen")
        ? static_cast<int16_t>(j["stage3_gapopen"].asInt())
        : INT16_MIN;
    sreq.stage3_gapext = j.isMember("stage3_gapext")
        ? static_cast<int16_t>(j["stage3_gapext"].asInt())
        : INT16_MIN;
    sreq.stage3_min_pident_x100 = static_cast<uint16_t>(j.get("stage3_min_pident_x100", 0).asUInt());
    sreq.stage3_min_nident = j.get("stage3_min_nident", 0).asUInt();
    sreq.context_abs = j.get("context_abs", 0).asUInt();
    sreq.context_frac_x10000 = static_cast<uint16_t>(j.get("context_frac_x10000", 0).asUInt());
    sreq.db_name = j.get("db_name", "").asString();

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

    // Phase 1: cached validation (synchronous, on Drogon thread)
    {
        std::lock_guard<std::mutex> lock(cache_mutex_);
        if (cache_valid_) {
            std::string err = validate_info(cached_info_, sreq.db_name,
                                            sreq.k, sreq.mode, false);
            if (!err.empty()) {
                callback(make_error_response(drogon::k400BadRequest, err));
                return;
            }
        }
    }

    // Offload blocking backend I/O to a worker thread.
    // Drogon event loop threads must not be blocked.
    auto backend = backend_;
    auto self = this;
    auto cb = std::make_shared<std::function<void(const drogon::HttpResponsePtr&)>>(
        std::move(callback));

    std::thread([self, backend, sreq = std::move(sreq), cb]() {
        // Phase 2: fresh validation with slot check
        {
            InfoResponse fresh;
            std::string info_err;
            if (backend->info(fresh, info_err)) {
                self->update_cache(fresh);
                std::string err = validate_info(fresh, sreq.db_name,
                                                sreq.k, sreq.mode, true);
                if (!err.empty()) {
                    (*cb)(make_error_response(drogon::k400BadRequest, err));
                    return;
                }
            }
        }

        // Phase 3: actual search
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
        result["stage1_score"] = sresp.stage1_score;
        if (sresp.stage3_traceback)
            result["stage3_traceback"] = sresp.stage3_traceback;
        const char* s1name = (sresp.stage1_score == 2) ? "matchscore" : "coverscore";

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
                }
                hobj["q_len"] = hit.q_length;
                if (sresp.mode != 1) {
                    hobj["s_start"] = hit.s_start;
                    hobj["s_end"] = hit.s_end;
                }
                hobj["s_len"] = hit.s_length;
                hobj[s1name] = hit.stage1_score;
                if (sresp.mode != 1) {
                    hobj["chainscore"] = hit.score;
                }
                if (sresp.mode == 3) {
                    hobj["alnscore"] = hit.alnscore;
                    if (sresp.stage3_traceback) {
                        hobj["pident"] = static_cast<double>(hit.pident_x100) / 100.0;
                        hobj["nident"] = hit.nident;
                        hobj["nmismatch"] = hit.nmismatch;
                        hobj["cigar"] = hit.cigar;
                        hobj["q_seq"] = hit.q_seq;
                        hobj["s_seq"] = hit.s_seq;
                    }
                }
                hobj["volume"] = hit.volume;
                hits_arr.append(std::move(hobj));
            }
            qobj["hits"] = std::move(hits_arr);
            if (qr.skipped != 0) {
                qobj["skipped"] = true;
            }
            if (qr.warnings != 0) {
                Json::Value warn_arr(Json::arrayValue);
                if (qr.warnings & kWarnMultiDegen) {
                    warn_arr.append("multi_degen");
                }
                qobj["warnings"] = std::move(warn_arr);
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
        result["max_active_sequences"] = iresp.max_active_sequences;
        result["active_sequences"] = iresp.active_sequences;

        Json::Value databases_arr(Json::arrayValue);
        for (const auto& db : iresp.databases) {
            Json::Value dbobj;
            dbobj["name"] = db.name;
            dbobj["default_k"] = db.default_k;
            dbobj["max_mode"] = db.max_mode;

            Json::Value groups_arr(Json::arrayValue);
            for (const auto& g : db.groups) {
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
            dbobj["kmer_groups"] = std::move(groups_arr);
            databases_arr.append(std::move(dbobj));
        }
        result["databases"] = std::move(databases_arr);

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
