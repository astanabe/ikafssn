#include "ikafssnhttpd/search_request_json.hpp"

#include <climits>
#include <memory>

#include <json/json.h>

namespace ikafssn {

bool parse_search_request_json(const std::string& body,
                               std::string& job_id,
                               SearchRequest& sreq,
                               std::string& error_msg) {
    Json::CharReaderBuilder reader_builder;
    std::unique_ptr<Json::CharReader> reader(reader_builder.newCharReader());
    Json::Value root;
    std::string parse_errors;
    if (!reader->parse(body.c_str(), body.c_str() + body.size(),
                       &root, &parse_errors)) {
        error_msg = "Invalid JSON body: " + parse_errors;
        return false;
    }
    if (!root.isObject()) {
        error_msg = "JSON body must be an object";
        return false;
    }

    if (!root.isMember("job_id") || !root["job_id"].isString()) {
        error_msg = "Missing or non-string 'job_id'";
        return false;
    }
    job_id = root["job_id"].asString();
    if (job_id.empty()) {
        error_msg = "'job_id' must not be empty";
        return false;
    }

    const auto& j = root;
    sreq = SearchRequest{};
    sreq.k = static_cast<uint8_t>(j.get("k", 0).asUInt());
    sreq.stage2_min_score = static_cast<uint16_t>(j.get("stage2_min_score", 0).asUInt());
    if (j.isMember("has_stage2_min_score") && j["has_stage2_min_score"].asBool()) {
        sreq.has_stage2_min_score = 1;
    }
    sreq.stage2_max_gap = static_cast<uint16_t>(j.get("stage2_max_gap", 0).asUInt());
    sreq.stage2_max_lookback = static_cast<uint16_t>(j.get("stage2_max_lookback", 0).asUInt());
    sreq.stage2_max_nhit_per_subject = static_cast<uint16_t>(j.get("stage2_max_nhit_per_subject", 0).asUInt());
    if (j.isMember("stage1_min_score_frac") && j["stage1_min_score_frac"].isDouble()) {
        double frac = j["stage1_min_score_frac"].asDouble();
        if (frac > 0 && frac < 1.0) {
            sreq.stage1_min_score_frac_x10000 = static_cast<uint16_t>(frac * 10000.0);
        }
    }
    sreq.stage2_min_nhit_diag = static_cast<uint8_t>(j.get("stage2_min_nhit_diag", 0).asUInt());
    sreq.stage1_topn = static_cast<uint16_t>(j.get("stage1_topn", 0).asUInt());
    sreq.stage1_min_score = static_cast<uint16_t>(j.get("stage1_min_score", 0).asUInt());
    sreq.nresult = static_cast<uint16_t>(j.get("nresult", 0).asUInt());
    sreq.mode = static_cast<uint8_t>(j.get("mode", 0).asUInt());
    sreq.stage1_score = static_cast<uint8_t>(j.get("stage1_score", 0).asUInt());
    sreq.accept_qdegen = static_cast<uint8_t>(j.get("accept_qdegen", 1).asUInt());
    sreq.strand = static_cast<int8_t>(j.get("strand", 0).asInt());

    sreq.stage3_traceback = static_cast<uint8_t>(j.get("stage3_traceback", 0).asUInt());
    sreq.stage3_gapopen = j.isMember("stage3_gapopen")
        ? static_cast<int16_t>(j["stage3_gapopen"].asInt()) : INT16_MIN;
    sreq.stage3_gapext = j.isMember("stage3_gapext")
        ? static_cast<int16_t>(j["stage3_gapext"].asInt()) : INT16_MIN;
    sreq.stage3_min_ppositive_x100 = static_cast<uint16_t>(j.get("stage3_min_ppositive_x100", 0).asUInt());
    sreq.stage3_min_npositive = j.get("stage3_min_npositive", 0).asUInt();
    sreq.context_abs = j.get("context_abs", 0).asUInt();
    sreq.context_frac_x10000 = static_cast<uint16_t>(j.get("context_frac_x10000", 0).asUInt());
    sreq.max_degen_expand = static_cast<uint16_t>(j.get("max_degen_expand", 0).asUInt());

    if (j.isMember("stage3_score_matrix")) {
        auto val = j["stage3_score_matrix"];
        if (val.isString()) {
            std::string sm = val.asString();
            if (sm == "degmatch") sreq.score_matrix = 1;
            else if (sm == "dnafull") sreq.score_matrix = 2;
            else if (sm == "nuc44") sreq.score_matrix = 3;
        } else if (val.isInt()) {
            sreq.score_matrix = static_cast<uint8_t>(val.asInt());
        }
    }
    if (j.isMember("t")) {
        sreq.t = static_cast<uint8_t>(j["t"].asInt());
    }
    if (j.isMember("template_type")) {
        auto val = j["template_type"];
        if (val.isString()) {
            std::string s = val.asString();
            if (s == "coding") sreq.template_type = 1;
            else if (s == "optimal") sreq.template_type = 2;
            else if (s == "both") sreq.template_type = 3;
        } else if (val.isInt()) {
            sreq.template_type = static_cast<uint8_t>(val.asInt());
        }
    }
    sreq.min_seq_length = j.get("min_seq_length", 0).asUInt();
    sreq.min_length_split = j.get("min_length_split", 0).asUInt();
    sreq.overlap_length = j.get("overlap_length", 0).asUInt();
    sreq.max_freq_build = j.get("max_freq_build", 1).asUInt64();
    sreq.db = j.get("db", "").asString();

    std::string mode_str = j.get("seqidlist_mode", "none").asString();
    if (mode_str == "include") sreq.seqidlist_mode = SeqidlistMode::kInclude;
    else if (mode_str == "exclude") sreq.seqidlist_mode = SeqidlistMode::kExclude;
    else sreq.seqidlist_mode = SeqidlistMode::kNone;

    if (j.isMember("seqids") && j["seqids"].isArray()) {
        for (const auto& s : j["seqids"]) {
            sreq.seqids.push_back(s.asString());
        }
    }

    if (!j.isMember("queries") || !j["queries"].isArray() ||
        j["queries"].empty()) {
        error_msg = "Missing or empty 'queries' array";
        return false;
    }
    for (const auto& q : j["queries"]) {
        if (!q.isMember("qseqid") || !q.isMember("sequence")) {
            error_msg = "Each query must have 'qseqid' and 'sequence'";
            return false;
        }
        QueryEntry entry;
        entry.qseqid = q["qseqid"].asString();
        entry.sequence = q["sequence"].asString();
        if (entry.sequence.empty()) {
            error_msg = "Query sequence must not be empty";
            return false;
        }
        sreq.queries.push_back(std::move(entry));
    }
    return true;
}

} // namespace ikafssn
