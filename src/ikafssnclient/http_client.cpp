#include "ikafssnclient/http_client.hpp"

#include <climits>
#include <curl/curl.h>
#include <json/json.h>

#include <memory>
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

std::string build_request_json(const SearchRequest& req) {
    Json::Value root;

    root["k"] = req.k;
    root["stage2_min_score"] = req.stage2_min_score;
    if (req.has_stage2_min_score) {
        root["has_stage2_min_score"] = true;
    }
    root["stage2_max_gap"] = req.stage2_max_gap;
    root["stage2_max_lookback"] = req.stage2_max_lookback;
    root["stage2_max_nhit_per_subject"] = req.stage2_max_nhit_per_subject;
    root["stage2_max_nhit_per_subject_mode"] = req.stage2_max_nhit_per_subject_mode;
    root["stage2_max_nhit_in_total"] = req.stage2_max_nhit_in_total;
    root["stage1_max_nhit_per_subject"] = req.stage1_max_nhit_per_subject;
    root["stage1_max_nhit_per_subject_mode"] = req.stage1_max_nhit_per_subject_mode;
    if (req.stage1_min_score_frac_x10000 != 0) {
        root["stage1_min_score_frac"] =
            static_cast<double>(req.stage1_min_score_frac_x10000) / 10000.0;
    }
    root["stage2_min_nhit_diag"] = req.stage2_min_nhit_diag;
    root["stage1_max_nhit_per_volume"] = req.stage1_max_nhit_per_volume;
    root["stage1_max_nhit_in_total"] = req.stage1_max_nhit_in_total;
    root["stage1_min_score"] = req.stage1_min_score;
    root["mode"] = req.mode;
    root["accept_qdegen"] = req.accept_qdegen;
    root["strand"] = req.strand;

    root["stage3_traceback"] = req.stage3_traceback;
    if (req.stage3_gapopen != INT16_MIN)
        root["stage3_gapopen"] = req.stage3_gapopen;
    if (req.stage3_gapext != INT16_MIN)
        root["stage3_gapext"] = req.stage3_gapext;
    root["stage3_min_ppositive_x100"] = req.stage3_min_ppositive_x100;
    root["stage3_min_npositive"] = req.stage3_min_npositive;
    root["stage3_max_nhit_per_subject"] = req.stage3_max_nhit_per_subject;
    root["stage3_max_nhit_per_subject_mode"] = req.stage3_max_nhit_per_subject_mode;
    root["stage3_max_nhit_in_total"] = req.stage3_max_nhit_in_total;
    root["context_abs"] = req.context_abs;
    root["context_frac_x10000"] = req.context_frac_x10000;
    root["max_degen_expand"] = req.max_degen_expand;
    if (req.score_matrix != 0) root["stage3_score_matrix"] = req.score_matrix;
    if (req.t > 0) root["t"] = req.t;
    if (req.template_type > 0) root["template_type"] = req.template_type;
    root["min_seq_length"] = req.min_seq_length;
    root["min_length_split"] = req.min_length_split;
    root["overlap_length"] = req.overlap_length;
    root["max_freq_build"] = static_cast<Json::UInt64>(req.max_freq_build);
    if (!req.db.empty()) root["db"] = req.db;

    switch (req.seqidlist_mode) {
    case SeqidlistMode::kInclude: root["seqidlist_mode"] = "include"; break;
    case SeqidlistMode::kExclude: root["seqidlist_mode"] = "exclude"; break;
    default:                       root["seqidlist_mode"] = "none";    break;
    }

    Json::Value seqids(Json::arrayValue);
    for (const auto& s : req.seqids) seqids.append(s);
    root["seqids"] = std::move(seqids);

    Json::Value queries(Json::arrayValue);
    for (const auto& q : req.queries) {
        Json::Value qobj;
        qobj["qseqid"] = q.qseqid;
        qobj["sequence"] = q.sequence;
        queries.append(std::move(qobj));
    }
    root["queries"] = std::move(queries);

    Json::StreamWriterBuilder writer;
    writer["indentation"] = "";
    return Json::writeString(writer, root);
}

bool parse_info_json(const std::string& body, InfoResponse& resp,
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

    if (root.isMember("error")) {
        error_msg = "Server error: " + root["error"].asString();
        return false;
    }

    resp.status = (root.get("status", "").asString() == "success") ? 0 : 1;
    resp.default_k = static_cast<uint8_t>(root.get("default_k", 0).asUInt());
    resp.max_queue_size = root.get("max_queue_size", 0).asInt();
    resp.queue_depth = root.get("queue_depth", 0).asInt();
    resp.max_nseq_per_req = root.get("max_nseq_per_req", 0).asInt();

    if (root.isMember("databases") && root["databases"].isArray()) {
        for (const auto& dbj : root["databases"]) {
            DatabaseInfo db;
            db.name = dbj.get("name", "").asString();
            db.default_k = static_cast<uint8_t>(dbj.get("default_k", 0).asUInt());
            db.max_mode = static_cast<uint8_t>(dbj.get("max_mode", 2).asUInt());

            if (dbj.isMember("kmer_groups") && dbj["kmer_groups"].isArray()) {
                for (const auto& gj : dbj["kmer_groups"]) {
                    KmerGroupInfo g;
                    g.k = static_cast<uint8_t>(gj.get("k", 0).asUInt());
                    std::string ktype = gj.get("kmer_type", "uint16").asString();
                    g.kmer_type = (ktype == "uint32") ? 1 : 0;
                    if (gj.isMember("t"))             g.t = static_cast<uint8_t>(gj["t"].asUInt());
                    if (gj.isMember("template_type")) g.template_type = static_cast<uint8_t>(gj["template_type"].asUInt());
                    g.min_seq_length = gj.get("min_seq_length", 0).asUInt();
                    g.min_length_split = gj.get("min_length_split", 0).asUInt();
                    g.overlap_length = gj.get("overlap_length", 0).asUInt();
                    g.max_freq_build = gj.get("max_freq_build", 1).asUInt64();
                    g.max_degen_expand = gj.get("max_degen_expand", 0).asUInt();

                    if (gj.isMember("volumes") && gj["volumes"].isArray()) {
                        for (const auto& vj : gj["volumes"]) {
                            VolumeInfo v;
                            v.volume_index = static_cast<uint16_t>(vj.get("volume_index", 0).asUInt());
                            v.num_sequences = vj.get("num_sequences", 0).asUInt();
                            v.total_postings = vj.get("total_postings", 0).asUInt64();
                            v.total_bases = vj.get("total_bases", 0).asUInt64();
                            v.db = vj.get("db", "").asString();
                            g.volumes.push_back(std::move(v));
                        }
                    }
                    db.groups.push_back(std::move(g));
                }
            }
            resp.databases.push_back(std::move(db));
        }
    }
    return true;
}

bool http_info(const std::string& base_url, InfoResponse& resp,
               std::string& error_msg, const HttpAuthConfig& auth) {
    std::string url = base_url;
    if (!url.empty() && url.back() == '/') url.pop_back();
    url += "/api/v1/info";

    CURL* curl = curl_easy_init();
    if (!curl) {
        error_msg = "Failed to initialize libcurl";
        return false;
    }

    std::string response_body;
    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_HTTPGET, 1L);
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_callback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response_body);

    if (!auth.userpwd.empty()) {
        curl_easy_setopt(curl, CURLOPT_HTTPAUTH, CURLAUTH_BASIC);
        curl_easy_setopt(curl, CURLOPT_USERPWD, auth.userpwd.c_str());
    }
    if (!auth.netrc_file.empty()) {
        curl_easy_setopt(curl, CURLOPT_NETRC, CURL_NETRC_OPTIONAL);
        curl_easy_setopt(curl, CURLOPT_NETRC_FILE, auth.netrc_file.c_str());
    }

    CURLcode res = curl_easy_perform(curl);
    if (res != CURLE_OK) {
        error_msg = "HTTP request failed: ";
        error_msg += curl_easy_strerror(res);
        curl_easy_cleanup(curl);
        return false;
    }
    long http_code = 0;
    curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &http_code);
    curl_easy_cleanup(curl);
    if (http_code != 200) {
        Json::CharReaderBuilder rb;
        std::unique_ptr<Json::CharReader> r(rb.newCharReader());
        Json::Value err_json;
        std::string errs;
        if (r->parse(response_body.c_str(),
                     response_body.c_str() + response_body.size(),
                     &err_json, &errs) && err_json.isMember("error")) {
            error_msg = "HTTP " + std::to_string(http_code) + ": " +
                        err_json["error"].asString();
        } else {
            error_msg = "HTTP " + std::to_string(http_code);
        }
        return false;
    }
    return parse_info_json(response_body, resp, error_msg);
}

} // namespace ikafssn
