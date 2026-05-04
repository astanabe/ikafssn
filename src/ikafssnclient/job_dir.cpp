#include "ikafssnclient/job_dir.hpp"
#include "ikafssnclient/zstd_io.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <random>
#include <sstream>

#include <json/json.h>

namespace fs = std::filesystem;

namespace ikafssn {

std::string default_jobs_root() { return ".ikafssnclient"; }

bool ensure_jobs_root(std::string& error_msg) {
    std::error_code ec;
    fs::create_directories(default_jobs_root(), ec);
    if (ec) {
        error_msg = "ensure_jobs_root: " + ec.message();
        return false;
    }
    return true;
}

static bool atomic_write_file(const std::string& path,
                              const std::string& body,
                              std::string& error_msg) {
    fs::path p(path);
    if (p.has_parent_path()) {
        std::error_code ec;
        fs::create_directories(p.parent_path(), ec);
    }
    std::string tmp = path + ".tmp";
    {
        std::ofstream out(tmp, std::ios::binary | std::ios::trunc);
        if (!out.is_open()) {
            error_msg = "atomic_write: cannot open " + tmp;
            return false;
        }
        out << body;
        if (!out.good()) {
            error_msg = "atomic_write: write failed for " + tmp;
            return false;
        }
    }
    if (std::rename(tmp.c_str(), path.c_str()) != 0) {
        error_msg = "atomic_write: rename failed for " + path;
        std::remove(tmp.c_str());
        return false;
    }
    return true;
}

static std::string group_dir(const std::string& root,
                             const std::string& group_id) {
    return root + "/" + group_id;
}

static std::string group_meta_path(const std::string& root,
                                   const std::string& group_id) {
    return group_dir(root, group_id) + "/group.json";
}

static std::string job_meta_path(const std::string& root,
                                 const std::string& group_id,
                                 const std::string& job_id) {
    return group_dir(root, group_id) + "/" + job_id + ".json";
}

static std::string job_deflines_path(const std::string& root,
                                     const std::string& group_id,
                                     const std::string& job_id) {
    return group_dir(root, group_id) + "/" + job_id + ".deflines.zst";
}

static std::string job_result_path(const std::string& root,
                                   const std::string& group_id,
                                   const std::string& job_id) {
    return group_dir(root, group_id) + "/" + job_id + ".result.bin.zst";
}

bool write_group_meta(const std::string& root, const GroupMeta& meta,
                      std::string& error_msg) {
    Json::Value root_j;
    root_j["group_id"]            = meta.group_id;
    root_j["submitted_at"]        = static_cast<Json::Int64>(meta.submitted_at);
    root_j["httpd_url"]           = meta.httpd_url;
    root_j["db"]                  = meta.db;
    root_j["query_file_path_abs"] = meta.query_file_path_abs;
    root_j["query_file_sha256"]   = meta.query_file_sha256;
    root_j["max_seqs_per_req"]    = meta.max_seqs_per_req;
    root_j["k"]                   = meta.k;
    root_j["mode"]                = meta.mode;
    root_j["t"]                   = meta.t;
    root_j["template_type"]       = meta.template_type;
    root_j["outfmt"]              = meta.outfmt;
    root_j["output_path"]         = meta.output_path;
    Json::Value jids(Json::arrayValue);
    for (const auto& j : meta.job_ids) jids.append(j);
    root_j["job_ids"] = std::move(jids);
    Json::Value sc;
    sc["queued"]  = meta.cnt_queued;
    sc["running"] = meta.cnt_running;
    sc["done"]    = meta.cnt_done;
    sc["failed"]  = meta.cnt_failed;
    root_j["status_counts"] = std::move(sc);

    Json::StreamWriterBuilder w;
    w["indentation"] = "  ";
    return atomic_write_file(group_meta_path(root, meta.group_id),
                             Json::writeString(w, root_j) + "\n",
                             error_msg);
}

static bool parse_json_file(const std::string& path, Json::Value& out,
                            std::string& error_msg) {
    std::ifstream in(path, std::ios::binary);
    if (!in.is_open()) {
        error_msg = "cannot open " + path;
        return false;
    }
    std::ostringstream oss;
    oss << in.rdbuf();
    std::string body = oss.str();
    Json::CharReaderBuilder rb;
    std::unique_ptr<Json::CharReader> r(rb.newCharReader());
    std::string errs;
    if (!r->parse(body.c_str(), body.c_str() + body.size(), &out, &errs)) {
        error_msg = "JSON parse failed for " + path + ": " + errs;
        return false;
    }
    return true;
}

bool read_group_meta(const std::string& root, const std::string& group_id,
                     GroupMeta& out, std::string& error_msg) {
    Json::Value j;
    if (!parse_json_file(group_meta_path(root, group_id), j, error_msg)) {
        return false;
    }
    out = GroupMeta{};
    out.group_id            = j.get("group_id", "").asString();
    out.submitted_at        = j.get("submitted_at", 0).asInt64();
    out.httpd_url           = j.get("httpd_url", "").asString();
    out.db                  = j.get("db", "").asString();
    out.query_file_path_abs = j.get("query_file_path_abs", "").asString();
    out.query_file_sha256   = j.get("query_file_sha256", "").asString();
    out.max_seqs_per_req    = j.get("max_seqs_per_req", 0).asInt();
    out.k                   = static_cast<uint8_t>(j.get("k", 0).asUInt());
    out.mode                = static_cast<uint8_t>(j.get("mode", 0).asUInt());
    out.t                   = static_cast<uint8_t>(j.get("t", 0).asUInt());
    out.template_type       = static_cast<uint8_t>(j.get("template_type", 0).asUInt());
    out.outfmt              = j.get("outfmt", "").asString();
    out.output_path         = j.get("output_path", "").asString();
    if (j.isMember("job_ids") && j["job_ids"].isArray()) {
        for (const auto& v : j["job_ids"]) out.job_ids.push_back(v.asString());
    }
    if (j.isMember("status_counts")) {
        const auto& sc = j["status_counts"];
        out.cnt_queued  = sc.get("queued", 0).asInt();
        out.cnt_running = sc.get("running", 0).asInt();
        out.cnt_done    = sc.get("done", 0).asInt();
        out.cnt_failed  = sc.get("failed", 0).asInt();
    }
    return true;
}

bool write_job_meta(const std::string& root, const JobMeta& meta,
                    std::string& error_msg) {
    Json::Value j;
    j["job_id"]         = meta.job_id;
    j["group_id"]       = meta.group_id;
    j["n_seqs"]         = meta.n_seqs;
    j["fasta_file"]     = meta.fasta_file;
    Json::Value range(Json::arrayValue);
    range.append(meta.seq_index_range.first);
    range.append(meta.seq_index_range.second);
    j["seq_index_range"] = std::move(range);
    j["status"]         = meta.status;
    j["attempts"]       = meta.attempts;
    j["submitted_at"]   = static_cast<Json::Int64>(meta.submitted_at);
    j["completed_at"]   = static_cast<Json::Int64>(meta.completed_at);
    j["last_polled_at"] = static_cast<Json::Int64>(meta.last_polled_at);
    if (!meta.error_message.empty()) j["error_message"] = meta.error_message;
    if (!meta.fail_reason.empty())   j["fail_reason"]   = meta.fail_reason;

    Json::StreamWriterBuilder w;
    w["indentation"] = "  ";
    return atomic_write_file(job_meta_path(root, meta.group_id, meta.job_id),
                             Json::writeString(w, j) + "\n", error_msg);
}

bool read_job_meta(const std::string& root, const std::string& group_id,
                   const std::string& job_id, JobMeta& out,
                   std::string& error_msg) {
    Json::Value j;
    if (!parse_json_file(job_meta_path(root, group_id, job_id),
                         j, error_msg)) return false;
    out = JobMeta{};
    out.job_id     = j.get("job_id", "").asString();
    out.group_id   = j.get("group_id", "").asString();
    out.n_seqs     = j.get("n_seqs", 0).asInt();
    out.fasta_file = j.get("fasta_file", "").asString();
    if (j.isMember("seq_index_range") && j["seq_index_range"].isArray()
        && j["seq_index_range"].size() == 2) {
        out.seq_index_range.first  = j["seq_index_range"][0].asInt();
        out.seq_index_range.second = j["seq_index_range"][1].asInt();
    }
    out.status         = j.get("status", "").asString();
    out.attempts       = j.get("attempts", 0).asInt();
    out.submitted_at   = j.get("submitted_at", 0).asInt64();
    out.completed_at   = j.get("completed_at", 0).asInt64();
    out.last_polled_at = j.get("last_polled_at", 0).asInt64();
    out.error_message  = j.get("error_message", "").asString();
    out.fail_reason    = j.get("fail_reason", "").asString();
    return true;
}

std::vector<std::string> list_groups(const std::string& root) {
    std::vector<std::pair<int64_t, std::string>> by_time;
    std::error_code ec;
    if (!fs::exists(root, ec)) return {};
    for (const auto& ent : fs::directory_iterator(root, ec)) {
        if (!ent.is_directory()) continue;
        std::string gid = ent.path().filename().string();
        GroupMeta gm;
        std::string err;
        if (!read_group_meta(root, gid, gm, err)) continue;
        by_time.emplace_back(gm.submitted_at, gid);
    }
    std::sort(by_time.begin(), by_time.end());
    std::vector<std::string> out;
    out.reserve(by_time.size());
    for (auto& [_, g] : by_time) out.push_back(std::move(g));
    return out;
}

std::vector<std::string> list_jobs(const std::string& root,
                                   const std::string& group_id) {
    GroupMeta gm;
    std::string err;
    if (!read_group_meta(root, group_id, gm, err)) return {};
    return gm.job_ids;
}

bool resolve_id(const std::string& root, const std::string& id,
                bool& is_group, std::string& group_id_out,
                std::string& job_id_out) {
    is_group = false;
    group_id_out.clear();
    job_id_out.clear();

    std::error_code ec;
    if (fs::exists(group_meta_path(root, id), ec)) {
        is_group = true;
        group_id_out = id;
        return true;
    }
    // Otherwise scan groups for a child with this job_id.
    auto groups = list_groups(root);
    for (const auto& g : groups) {
        if (fs::exists(job_meta_path(root, g, id), ec)) {
            group_id_out = g;
            job_id_out = id;
            return true;
        }
    }
    return false;
}

bool write_job_deflines(const std::string& root,
                        const std::string& group_id,
                        const std::string& job_id,
                        const std::vector<std::string>& deflines,
                        std::string& error_msg) {
    std::string body;
    for (const auto& s : deflines) {
        body.append(s);
        body.push_back('\n');
    }
    fs::path p = job_deflines_path(root, group_id, job_id);
    std::error_code ec;
    fs::create_directories(p.parent_path(), ec);
    return zstd_compress_to_file(p.string(), body.data(), body.size(),
                                  3, error_msg);
}

bool read_job_deflines(const std::string& root,
                       const std::string& group_id,
                       const std::string& job_id,
                       std::vector<std::string>& out,
                       std::string& error_msg) {
    std::vector<uint8_t> blob;
    if (!zstd_decompress_file(job_deflines_path(root, group_id, job_id),
                              blob, error_msg)) return false;
    std::string body(blob.begin(), blob.end());
    std::string line;
    for (char c : body) {
        if (c == '\n') { out.push_back(line); line.clear(); }
        else            { line.push_back(c); }
    }
    if (!line.empty()) out.push_back(std::move(line));
    return true;
}

bool write_job_result(const std::string& root,
                      const std::string& group_id,
                      const std::string& job_id,
                      const std::vector<uint8_t>& blob,
                      std::string& error_msg) {
    fs::path p = job_result_path(root, group_id, job_id);
    std::error_code ec;
    fs::create_directories(p.parent_path(), ec);
    return zstd_compress_to_file(p.string(), blob.data(), blob.size(),
                                  3, error_msg);
}

bool read_job_result(const std::string& root,
                     const std::string& group_id,
                     const std::string& job_id,
                     std::vector<uint8_t>& out,
                     std::string& error_msg) {
    return zstd_decompress_file(job_result_path(root, group_id, job_id),
                                out, error_msg);
}

bool delete_job_result(const std::string& root,
                       const std::string& group_id,
                       const std::string& job_id,
                       std::string& error_msg) {
    std::error_code ec;
    fs::remove(job_result_path(root, group_id, job_id), ec);
    if (ec) {
        error_msg = "delete_job_result: " + ec.message();
        return false;
    }
    return true;
}

std::string make_uuidv4() {
    static thread_local std::mt19937_64 rng(
        std::chrono::steady_clock::now().time_since_epoch().count() ^
        static_cast<uint64_t>(reinterpret_cast<uintptr_t>(&rng)));
    uint64_t a = rng();
    uint64_t b = rng();
    // Set version (4) and variant (10xx).
    uint16_t time_hi = static_cast<uint16_t>((a >> 16) & 0x0FFFu) | 0x4000u;
    uint16_t clock_seq = static_cast<uint16_t>((b >> 48) & 0x3FFFu) | 0x8000u;
    char buf[40];
    std::snprintf(buf, sizeof(buf),
                  "%08x-%04x-%04x-%04x-%012llx",
                  static_cast<unsigned>(a >> 32),
                  static_cast<unsigned>(a & 0xFFFFu),
                  static_cast<unsigned>(time_hi),
                  static_cast<unsigned>(clock_seq),
                  static_cast<unsigned long long>(b & 0xFFFFFFFFFFFFull));
    return buf;
}

} // namespace ikafssn
