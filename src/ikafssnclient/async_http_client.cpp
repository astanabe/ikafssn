#include "ikafssnclient/async_http_client.hpp"

#include <curl/curl.h>
#include <json/json.h>

#include <memory>
#include <sstream>

namespace ikafssn {

namespace {

size_t curl_write_string(char* ptr, size_t s, size_t n, void* ud) {
    auto* str = static_cast<std::string*>(ud);
    str->append(ptr, s * n);
    return s * n;
}

void apply_auth(CURL* c, const HttpAuthConfig& auth) {
    if (!auth.userpwd.empty()) {
        curl_easy_setopt(c, CURLOPT_HTTPAUTH, CURLAUTH_BASIC);
        curl_easy_setopt(c, CURLOPT_USERPWD, auth.userpwd.c_str());
    }
    if (!auth.netrc_file.empty()) {
        curl_easy_setopt(c, CURLOPT_NETRC, CURL_NETRC_OPTIONAL);
        curl_easy_setopt(c, CURLOPT_NETRC_FILE, auth.netrc_file.c_str());
    }
}

AsyncHttpOutcome classify(long http_code) {
    if (http_code >= 200 && http_code < 300) return AsyncHttpOutcome::kOk;
    if (http_code == 404) return AsyncHttpOutcome::kNotFound;
    if (http_code == 409) return AsyncHttpOutcome::kConflict;
    if (http_code == 400) return AsyncHttpOutcome::kBadRequest;
    if (http_code >= 500) return AsyncHttpOutcome::kServerError;
    return AsyncHttpOutcome::kTransport;
}

} // namespace

AsyncHttpOutcome http_submit_job(const std::string& base_url,
                                 const std::string& job_id,
                                 const SearchRequest& req,
                                 std::string& error_msg,
                                 const HttpAuthConfig& auth) {
    std::string url = base_url;
    if (!url.empty() && url.back() == '/') url.pop_back();
    url += "/api/v1/jobs";

    // Inject job_id into the JSON body produced by build_request_json.
    std::string base_body = build_request_json(req);
    Json::CharReaderBuilder rb;
    std::unique_ptr<Json::CharReader> reader(rb.newCharReader());
    Json::Value root;
    std::string parse_err;
    if (!reader->parse(base_body.c_str(), base_body.c_str() + base_body.size(),
                       &root, &parse_err)) {
        error_msg = "submit_job: built body could not be re-parsed: " + parse_err;
        return AsyncHttpOutcome::kTransport;
    }
    root["job_id"] = job_id;
    Json::StreamWriterBuilder writer;
    writer["indentation"] = "";
    std::string body = Json::writeString(writer, root);

    CURL* c = curl_easy_init();
    if (!c) {
        error_msg = "Failed to initialize libcurl";
        return AsyncHttpOutcome::kTransport;
    }

    std::string resp;
    curl_easy_setopt(c, CURLOPT_URL, url.c_str());
    curl_easy_setopt(c, CURLOPT_POST, 1L);
    curl_easy_setopt(c, CURLOPT_POSTFIELDS, body.c_str());
    curl_easy_setopt(c, CURLOPT_POSTFIELDSIZE, static_cast<long>(body.size()));
    curl_easy_setopt(c, CURLOPT_WRITEFUNCTION, curl_write_string);
    curl_easy_setopt(c, CURLOPT_WRITEDATA, &resp);

    curl_slist* hdrs = nullptr;
    hdrs = curl_slist_append(hdrs, "Content-Type: application/json");
    curl_easy_setopt(c, CURLOPT_HTTPHEADER, hdrs);
    apply_auth(c, auth);

    CURLcode rc = curl_easy_perform(c);
    long code = 0;
    curl_easy_getinfo(c, CURLINFO_RESPONSE_CODE, &code);
    curl_slist_free_all(hdrs);
    curl_easy_cleanup(c);

    if (rc != CURLE_OK) {
        error_msg = std::string("submit_job HTTP error: ")
                  + curl_easy_strerror(rc);
        return AsyncHttpOutcome::kTransport;
    }

    auto outcome = classify(code);
    if (outcome != AsyncHttpOutcome::kOk &&
        outcome != AsyncHttpOutcome::kConflict) {
        Json::Value err_json;
        std::string errs;
        if (reader->parse(resp.c_str(), resp.c_str() + resp.size(),
                          &err_json, &errs) && err_json.isMember("error")) {
            error_msg = "HTTP " + std::to_string(code) + ": " +
                        err_json["error"].asString();
        } else {
            error_msg = "HTTP " + std::to_string(code);
        }
    }
    return outcome;
}

AsyncHttpOutcome http_get_job_status(const std::string& base_url,
                                     const std::string& job_id,
                                     AsyncJobStatus& out,
                                     std::string& error_msg,
                                     const HttpAuthConfig& auth) {
    std::string url = base_url;
    if (!url.empty() && url.back() == '/') url.pop_back();
    url += "/api/v1/jobs/" + job_id;

    CURL* c = curl_easy_init();
    if (!c) {
        error_msg = "curl init";
        return AsyncHttpOutcome::kTransport;
    }
    std::string resp;
    curl_easy_setopt(c, CURLOPT_URL, url.c_str());
    curl_easy_setopt(c, CURLOPT_HTTPGET, 1L);
    curl_easy_setopt(c, CURLOPT_WRITEFUNCTION, curl_write_string);
    curl_easy_setopt(c, CURLOPT_WRITEDATA, &resp);
    apply_auth(c, auth);

    CURLcode rc = curl_easy_perform(c);
    long code = 0;
    curl_easy_getinfo(c, CURLINFO_RESPONSE_CODE, &code);
    curl_easy_cleanup(c);

    if (rc != CURLE_OK) {
        error_msg = std::string("get_status: ") + curl_easy_strerror(rc);
        return AsyncHttpOutcome::kTransport;
    }

    auto outcome = classify(code);
    if (outcome != AsyncHttpOutcome::kOk) {
        error_msg = "HTTP " + std::to_string(code);
        return outcome;
    }

    Json::CharReaderBuilder rb;
    std::unique_ptr<Json::CharReader> reader(rb.newCharReader());
    Json::Value root;
    std::string errs;
    if (!reader->parse(resp.c_str(), resp.c_str() + resp.size(),
                       &root, &errs)) {
        error_msg = "get_status: invalid JSON reply";
        return AsyncHttpOutcome::kTransport;
    }
    out = AsyncJobStatus{};
    out.status        = root.get("status", "").asString();
    out.error_message = root.get("error_message", "").asString();
    out.fail_reason   = root.get("fail_reason", "").asString();
    out.attempts      = root.get("attempts", 0).asInt();
    out.submitted_at  = root.get("submitted_at", 0).asInt64();
    out.started_at    = root.get("started_at", 0).asInt64();
    out.completed_at  = root.get("completed_at", 0).asInt64();
    return AsyncHttpOutcome::kOk;
}

AsyncHttpOutcome http_get_job_result(const std::string& base_url,
                                     const std::string& job_id,
                                     std::vector<uint8_t>& out,
                                     std::string& error_msg,
                                     const HttpAuthConfig& auth) {
    std::string url = base_url;
    if (!url.empty() && url.back() == '/') url.pop_back();
    url += "/api/v1/jobs/" + job_id + "/result";

    CURL* c = curl_easy_init();
    if (!c) {
        error_msg = "curl init";
        return AsyncHttpOutcome::kTransport;
    }
    std::string resp;
    curl_easy_setopt(c, CURLOPT_URL, url.c_str());
    curl_easy_setopt(c, CURLOPT_HTTPGET, 1L);
    curl_easy_setopt(c, CURLOPT_WRITEFUNCTION, curl_write_string);
    curl_easy_setopt(c, CURLOPT_WRITEDATA, &resp);
    apply_auth(c, auth);

    CURLcode rc = curl_easy_perform(c);
    long code = 0;
    curl_easy_getinfo(c, CURLINFO_RESPONSE_CODE, &code);
    curl_easy_cleanup(c);

    if (rc != CURLE_OK) {
        error_msg = std::string("get_result: ") + curl_easy_strerror(rc);
        return AsyncHttpOutcome::kTransport;
    }

    auto outcome = classify(code);
    if (outcome != AsyncHttpOutcome::kOk) {
        error_msg = "HTTP " + std::to_string(code);
        return outcome;
    }
    out.assign(resp.begin(), resp.end());
    return AsyncHttpOutcome::kOk;
}

} // namespace ikafssn
