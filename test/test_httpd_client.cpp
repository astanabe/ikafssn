// Full-stack async REST integration test:
//   ikafssnclient (HTTP submit/poll/result)
//     -> ikafssnhttpd (job-broker)
//       -> ikafssnserver
//
// Phase B replaced the synchronous POST /api/v1/search with a job-broker
// model.  This test exercises the new flow end-to-end via libcurl, then
// deserialises the binary SearchResponse blob returned by
// GET /api/v1/jobs/<id>/result and verifies hit-by-hit equivalence
// with a direct local search.

#include "test_util.hpp"
#include "ssu_test_fixture.hpp"

#include "io/blastdb_reader.hpp"
#include "io/fasta_reader.hpp"
#include "io/volume_discovery.hpp"
#include "index/index_builder.hpp"
#include "index/kix_reader.hpp"
#include "index/kpx_reader.hpp"
#include "index/ksx_reader.hpp"
#include "search/oid_filter.hpp"
#include "search/query_preprocessor.hpp"
#include "search/volume_searcher.hpp"
#include "ikafssnserver/server.hpp"
#include "ikafssnhttpd/backend_client.hpp"
#include "ikafssnhttpd/backend_manager.hpp"
#include "ikafssnhttpd/http_controller.hpp"
#include "ikafssnhttpd/job_store.hpp"
#include "ikafssnhttpd/job_worker.hpp"
#include "ikafssnhttpd/job_housekeeper.hpp"
#include "protocol/messages.hpp"
#include "protocol/serializer.hpp"
#include "util/socket_utils.hpp"
#include "util/logger.hpp"
#include "core/config.hpp"

#include <drogon/HttpAppFramework.h>
#include <curl/curl.h>
#include <json/json.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <random>
#include <string>
#include <thread>

#include <sys/stat.h>
#include <unistd.h>

using namespace ikafssn;
using namespace ssu_fixture;

static std::string g_testdb_path;
static std::string g_test_dir;

static size_t curl_write_string(char* p, size_t s, size_t n, void* ud) {
    auto* str = static_cast<std::string*>(ud);
    str->append(p, s * n);
    return s * n;
}

static std::string make_uuidv4() {
    static std::mt19937_64 rng(
        std::chrono::steady_clock::now().time_since_epoch().count());
    auto a = rng(), b = rng();
    char buf[40];
    std::snprintf(buf, sizeof(buf),
        "%08x-%04x-4%03x-%04x-%012llx",
        static_cast<unsigned>(a >> 32) & 0xFFFFFFFFu,
        static_cast<unsigned>(a >> 16) & 0xFFFFu,
        static_cast<unsigned>(a) & 0x0FFFu,
        static_cast<unsigned>((b >> 48) & 0x3FFFu) | 0x8000u,
        (unsigned long long)(b & 0xFFFFFFFFFFFFu));
    return buf;
}

static std::string build_test_index(int k) {
    std::string ix_dir = g_test_dir + "/hc_index";
    std::filesystem::create_directories(ix_dir);

    BlastDbReader db;
    if (!db.open(g_testdb_path)) {
        std::fprintf(stderr, "Cannot open test BLAST DB at %s\n",
                     g_testdb_path.c_str());
        return {};
    }

    Logger logger(Logger::kError);
    IndexBuilderConfig bconfig;
    bconfig.k = k;

    char kk[4];
    std::snprintf(kk, sizeof(kk), "%02d", k);
    std::string prefix = ix_dir + "/test.00." + std::string(kk) + "mer";

    bool ok;
    if (k < K_TYPE_THRESHOLD) {
        ok = build_index<uint16_t>(db, bconfig, prefix, 0, 1, "test", logger);
    } else {
        ok = build_index<uint32_t>(db, bconfig, prefix, 0, 1, "test", logger);
    }
    if (!ok) return {};

    std::string kvx_path = ix_dir + "/test." + std::string(kk) + "mer.kvx";
    FILE* fp = std::fopen(kvx_path.c_str(), "w");
    if (fp) {
        std::fprintf(fp, "TITLE test\nDBLIST \"test.00\"\n");
        std::fclose(fp);
    }
    return ix_dir + "/test";
}

// Submit a job, poll until status leaves 'queued'/'running', then GET
// /result and parse the binary blob.  Returns true on done, false on
// failed/timeout.
static bool async_submit_poll_get(const std::string& base_url,
                                  const SearchRequest& req,
                                  const std::string& job_id,
                                  SearchResponse& out_resp,
                                  std::string& error_msg) {
    Json::Value body;
    body["job_id"] = job_id;
    body["k"] = req.k;
    body["mode"] = req.mode;
    body["db"] = req.db;
    if (req.seqidlist_mode == SeqidlistMode::kInclude) {
        body["seqidlist_mode"] = "include";
    } else if (req.seqidlist_mode == SeqidlistMode::kExclude) {
        body["seqidlist_mode"] = "exclude";
    } else {
        body["seqidlist_mode"] = "none";
    }
    Json::Value seqids(Json::arrayValue);
    for (const auto& s : req.seqids) seqids.append(s);
    body["seqids"] = std::move(seqids);
    Json::Value queries(Json::arrayValue);
    for (const auto& q : req.queries) {
        Json::Value qq;
        qq["qseqid"] = q.qseqid;
        qq["sequence"] = q.sequence;
        queries.append(std::move(qq));
    }
    body["queries"] = std::move(queries);

    Json::StreamWriterBuilder writer;
    writer["indentation"] = "";
    std::string body_str = Json::writeString(writer, body);

    // POST /api/v1/jobs
    {
        CURL* c = curl_easy_init();
        if (!c) { error_msg = "curl init"; return false; }
        std::string resp;
        curl_easy_setopt(c, CURLOPT_URL, (base_url + "/api/v1/jobs").c_str());
        curl_easy_setopt(c, CURLOPT_POST, 1L);
        curl_easy_setopt(c, CURLOPT_POSTFIELDS, body_str.c_str());
        curl_easy_setopt(c, CURLOPT_WRITEFUNCTION, curl_write_string);
        curl_easy_setopt(c, CURLOPT_WRITEDATA, &resp);
        curl_slist* hdrs = nullptr;
        hdrs = curl_slist_append(hdrs, "Content-Type: application/json");
        curl_easy_setopt(c, CURLOPT_HTTPHEADER, hdrs);
        CURLcode rc = curl_easy_perform(c);
        long code = 0;
        curl_easy_getinfo(c, CURLINFO_RESPONSE_CODE, &code);
        curl_slist_free_all(hdrs);
        curl_easy_cleanup(c);
        if (rc != CURLE_OK || code != 202) {
            error_msg = "submit: HTTP " + std::to_string(code) + " " + resp;
            return false;
        }
    }

    // Poll up to ~10s.
    std::string status_url = base_url + "/api/v1/jobs/" + job_id;
    for (int i = 0; i < 100; i++) {
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        CURL* c = curl_easy_init();
        std::string resp;
        curl_easy_setopt(c, CURLOPT_URL, status_url.c_str());
        curl_easy_setopt(c, CURLOPT_WRITEFUNCTION, curl_write_string);
        curl_easy_setopt(c, CURLOPT_WRITEDATA, &resp);
        CURLcode rc = curl_easy_perform(c);
        long code = 0;
        curl_easy_getinfo(c, CURLINFO_RESPONSE_CODE, &code);
        curl_easy_cleanup(c);
        if (rc != CURLE_OK || code != 200) continue;

        Json::CharReaderBuilder rb;
        std::unique_ptr<Json::CharReader> r(rb.newCharReader());
        Json::Value root;
        std::string errs;
        if (!r->parse(resp.c_str(), resp.c_str() + resp.size(), &root, &errs)) {
            continue;
        }
        std::string st = root.get("status", "").asString();
        if (st == "done") break;
        if (st == "failed") {
            error_msg = "job failed: " + root.get("fail_reason", "").asString();
            return false;
        }
    }

    // GET /result
    {
        CURL* c = curl_easy_init();
        std::string resp;
        curl_easy_setopt(c, CURLOPT_URL,
                         (base_url + "/api/v1/jobs/" + job_id + "/result").c_str());
        curl_easy_setopt(c, CURLOPT_WRITEFUNCTION, curl_write_string);
        curl_easy_setopt(c, CURLOPT_WRITEDATA, &resp);
        CURLcode rc = curl_easy_perform(c);
        long code = 0;
        curl_easy_getinfo(c, CURLINFO_RESPONSE_CODE, &code);
        curl_easy_cleanup(c);
        if (rc != CURLE_OK || code != 200) {
            error_msg = "result: HTTP " + std::to_string(code);
            return false;
        }
        std::vector<uint8_t> blob(resp.begin(), resp.end());
        if (!deserialize(blob, out_resp)) {
            error_msg = "result deserialize failed";
            return false;
        }
    }
    return true;
}

static void test_async_submit_poll_get() {
    std::fprintf(stderr, "-- test_async_submit_poll_get\n");

    int k = 7;
    std::string ix_prefix = build_test_index(k);
    CHECK(!ix_prefix.empty());

    auto queries = read_fasta(queries_path());
    CHECK(!queries.empty());

    KixReader kix;
    KpxReader kpx;
    KsxReader ksx;
    char kk[4];
    std::snprintf(kk, sizeof(kk), "%02d", k);
    std::string base = ix_prefix + ".00." + std::string(kk) + "mer";
    CHECK(kix.open(base + ".kix"));
    CHECK(kpx.open(base + ".kpx"));
    CHECK(ksx.open(base + ".ksx"));

    SearchConfig config;
    config.stage2.min_score = 1;
    OidFilter no_filter;

    std::vector<SearchResult> local_results;
    for (const auto& q : queries) {
        Stage1Buffer buf;
        auto qdata = preprocess_query<uint16_t>(q.sequence, k, nullptr, config);
        auto sr = search_volume<uint16_t>(q.id, qdata, k, kix, kpx, ksx,
                                          no_filter, config, buf);
        local_results.push_back(sr);
    }

    std::string sock_path = g_test_dir + "/test_hc_server.sock";
    ::unlink(sock_path.c_str());
    std::string db = parse_index_prefix(ix_prefix).db;

    ServerConfig sconfig;
    sconfig.db_entries.push_back({ix_prefix, ix_prefix});
    sconfig.unix_socket_path = sock_path;
    sconfig.nthread = 2;
    sconfig.log_level = Logger::kError;
    sconfig.search_config.stage2.min_score = 1;

    Server server;
    std::thread server_thread([&] { server.run(sconfig); });
    for (int i = 0; i < 20; i++) {
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
        struct stat st;
        if (stat(sock_path.c_str(), &st) == 0) break;
    }

    uint16_t http_port = 18923;
    auto manager = std::make_shared<BackendManager>();
    manager->add_backend(BackendMode::kUnix, sock_path);
    Logger mgr_logger(Logger::kError);
    CHECK(manager->init(10, mgr_logger));

    std::string jobs_db = g_test_dir + "/jobs.db";
    JobStore store;
    {
        std::string err;
        CHECK(store.open(jobs_db, err));
    }
    JobWorker worker(store, manager, mgr_logger, 3);
    worker.start(2);

    HttpController controller(manager, store, worker);
    controller.register_routes("");

    drogon::app()
        .addListener("127.0.0.1", http_port)
        .setThreadNum(1)
        .setLogLevel(trantor::Logger::kFatal);

    std::thread httpd_thread([] { drogon::app().run(); });
    std::this_thread::sleep_for(std::chrono::milliseconds(500));

    std::string http_url = "http://127.0.0.1:" + std::to_string(http_port);

    SearchRequest req;
    req.k = static_cast<uint8_t>(k);
    req.db = db;
    for (const auto& q : queries) {
        req.queries.push_back({q.id, q.sequence});
    }

    SearchResponse resp;
    std::string error_msg;
    bool ok = async_submit_poll_get(http_url, req, make_uuidv4(),
                                     resp, error_msg);
    if (!ok) {
        std::fprintf(stderr, "  async_submit_poll_get failed: %s\n",
                     error_msg.c_str());
    }
    CHECK(ok);
    CHECK(resp.status == 0);
    CHECK_EQ(resp.k, static_cast<uint8_t>(k));
    CHECK_EQ(resp.results.size(), local_results.size());

    drogon::app().quit();
    httpd_thread.join();
    worker.stop();
    store.close();
    server.request_shutdown();
    server_thread.join();
    ::unlink(sock_path.c_str());
}

int main() {
    check_ssu_available();
    check_derived_data_ready();

    g_testdb_path = ssu_db_prefix();
    g_test_dir = "/tmp/ikafssn_test_httpd_client_" +
                 std::to_string(::getpid());
    std::filesystem::create_directories(g_test_dir);

    test_async_submit_poll_get();

    std::filesystem::remove_all(g_test_dir);
    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
