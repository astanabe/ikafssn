#include "ikafssnclient/http_client.hpp"

#include <curl/curl.h>
#include <json/json.h>

#include <sstream>

namespace ikafssn {

// libcurl write callback: append received data to a std::string.
static size_t write_callback(char* ptr, size_t size, size_t nmemb,
                             void* userdata) {
    auto* buf = static_cast<std::string*>(userdata);
    size_t total = size * nmemb;
    buf->append(ptr, total);
    return total;
}

// Build JSON request body from SearchRequest.
static std::string build_request_json(const SearchRequest& req) {
    Json::Value root;

    root["k"] = req.k;
    root["min_score"] = req.min_score;
    root["max_gap"] = req.max_gap;
    root["max_freq"] = req.max_freq;
    root["min_diag_hits"] = req.min_diag_hits;
    root["stage1_topn"] = req.stage1_topn;
    root["min_stage1_score"] = req.min_stage1_score;
    root["num_results"] = req.num_results;

    switch (req.seqidlist_mode) {
    case SeqidlistMode::kInclude:
        root["seqidlist_mode"] = "include";
        break;
    case SeqidlistMode::kExclude:
        root["seqidlist_mode"] = "exclude";
        break;
    default:
        root["seqidlist_mode"] = "none";
        break;
    }

    Json::Value seqids(Json::arrayValue);
    for (const auto& s : req.seqids) {
        seqids.append(s);
    }
    root["seqids"] = std::move(seqids);

    Json::Value queries(Json::arrayValue);
    for (const auto& q : req.queries) {
        Json::Value qobj;
        qobj["query_id"] = q.query_id;
        qobj["sequence"] = q.sequence;
        queries.append(std::move(qobj));
    }
    root["queries"] = std::move(queries);

    Json::StreamWriterBuilder writer;
    writer["indentation"] = "";
    return Json::writeString(writer, root);
}

// Parse JSON response into SearchResponse.
static bool parse_response_json(const std::string& body,
                                SearchResponse& resp,
                                std::string& error_msg) {
    Json::CharReaderBuilder reader_builder;
    std::unique_ptr<Json::CharReader> reader(reader_builder.newCharReader());

    Json::Value root;
    std::string parse_errors;
    if (!reader->parse(body.c_str(), body.c_str() + body.size(),
                       &root, &parse_errors)) {
        error_msg = "Failed to parse JSON response: " + parse_errors;
        return false;
    }

    // Check for error field
    if (root.isMember("error")) {
        error_msg = "Server error: " + root["error"].asString();
        return false;
    }

    std::string status = root.get("status", "").asString();
    resp.status = (status == "success") ? 0 : 1;
    resp.k = static_cast<uint8_t>(root.get("k", 0).asUInt());

    const auto& results = root["results"];
    if (!results.isArray()) {
        error_msg = "Missing 'results' array in response";
        return false;
    }

    for (const auto& qr : results) {
        QueryResult query_result;
        query_result.query_id = qr.get("query_id", "").asString();

        const auto& hits = qr["hits"];
        if (hits.isArray()) {
            for (const auto& h : hits) {
                ResponseHit hit;
                hit.accession = h.get("accession", "").asString();
                std::string strand = h.get("strand", "+").asString();
                hit.strand = (strand == "-") ? 1 : 0;
                hit.q_start = h.get("q_start", 0).asUInt();
                hit.q_end = h.get("q_end", 0).asUInt();
                hit.s_start = h.get("s_start", 0).asUInt();
                hit.s_end = h.get("s_end", 0).asUInt();
                hit.score = static_cast<uint16_t>(h.get("score", 0).asUInt());
                hit.volume = static_cast<uint16_t>(h.get("volume", 0).asUInt());
                query_result.hits.push_back(std::move(hit));
            }
        }

        resp.results.push_back(std::move(query_result));
    }

    return true;
}

bool http_search(const std::string& base_url, const SearchRequest& req,
                 SearchResponse& resp, std::string& error_msg) {
    // Build URL
    std::string url = base_url;
    if (!url.empty() && url.back() == '/') {
        url.pop_back();
    }
    url += "/api/v1/search";

    // Build JSON body
    std::string body = build_request_json(req);

    // Initialize curl
    CURL* curl = curl_easy_init();
    if (!curl) {
        error_msg = "Failed to initialize libcurl";
        return false;
    }

    std::string response_body;

    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_POST, 1L);
    curl_easy_setopt(curl, CURLOPT_POSTFIELDS, body.c_str());
    curl_easy_setopt(curl, CURLOPT_POSTFIELDSIZE,
                     static_cast<long>(body.size()));
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_callback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response_body);

    struct curl_slist* headers = nullptr;
    headers = curl_slist_append(headers, "Content-Type: application/json");
    curl_easy_setopt(curl, CURLOPT_HTTPHEADER, headers);

    // Perform request
    CURLcode res = curl_easy_perform(curl);

    if (res != CURLE_OK) {
        error_msg = "HTTP request failed: ";
        error_msg += curl_easy_strerror(res);
        curl_slist_free_all(headers);
        curl_easy_cleanup(curl);
        return false;
    }

    // Check HTTP status code
    long http_code = 0;
    curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &http_code);

    curl_slist_free_all(headers);
    curl_easy_cleanup(curl);

    if (http_code != 200) {
        // Try to parse error from JSON body
        Json::CharReaderBuilder rb;
        std::unique_ptr<Json::CharReader> r(rb.newCharReader());
        Json::Value err_json;
        std::string errs;
        if (r->parse(response_body.c_str(),
                     response_body.c_str() + response_body.size(),
                     &err_json, &errs) &&
            err_json.isMember("error")) {
            error_msg = "HTTP " + std::to_string(http_code) + ": " +
                        err_json["error"].asString();
        } else {
            error_msg = "HTTP " + std::to_string(http_code);
        }
        return false;
    }

    return parse_response_json(response_body, resp, error_msg);
}

} // namespace ikafssn
