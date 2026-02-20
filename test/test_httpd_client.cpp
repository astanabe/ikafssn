#include "test_util.hpp"

#include "io/blastdb_reader.hpp"
#include "io/fasta_reader.hpp"
#include "index/index_builder.hpp"
#include "index/kix_reader.hpp"
#include "index/kpx_reader.hpp"
#include "index/ksx_reader.hpp"
#include "search/oid_filter.hpp"
#include "search/volume_searcher.hpp"
#include "ikafssnserver/server.hpp"
#include "ikafssnhttpd/backend_client.hpp"
#include "ikafssnhttpd/http_controller.hpp"
#include "ikafssnclient/http_client.hpp"
#include "protocol/messages.hpp"
#include "util/socket_utils.hpp"
#include "util/logger.hpp"
#include "core/config.hpp"

#include <drogon/HttpAppFramework.h>

#include <cstdio>
#include <filesystem>
#include <string>
#include <thread>
#include <chrono>
#include <atomic>

#include <sys/stat.h>
#include <unistd.h>

using namespace ikafssn;

static std::string g_testdb_path;
static std::string g_test_dir;

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
    bconfig.partitions = 1;
    bconfig.buffer_size = uint64_t(1) << 30;

    char kk[4];
    std::snprintf(kk, sizeof(kk), "%02d", k);
    std::string prefix =
        ix_dir + "/test.00." + std::string(kk) + "mer";

    bool ok;
    if (k < K_TYPE_THRESHOLD) {
        ok = build_index<uint16_t>(db, bconfig, prefix, 0, 1, "test", logger);
    } else {
        ok = build_index<uint32_t>(db, bconfig, prefix, 0, 1, "test", logger);
    }

    if (!ok) {
        std::fprintf(stderr, "Failed to build index for k=%d\n", k);
        return {};
    }

    return ix_dir;
}

// Full stack test: ikafssnclient (http_search) -> ikafssnhttpd -> ikafssnserver
static void test_http_client_search() {
    std::fprintf(stderr, "-- test_http_client_search\n");

    int k = 7;
    std::string ix_dir = build_test_index(k);
    CHECK(!ix_dir.empty());

    // Read query FASTA
    std::string query_fasta =
        std::string(SOURCE_DIR) + "/test/testdata/queries.fasta";
    auto queries = read_fasta(query_fasta);
    CHECK(!queries.empty());

    // --- Direct local search for reference ---
    KixReader kix;
    KpxReader kpx;
    KsxReader ksx;
    char kk[4];
    std::snprintf(kk, sizeof(kk), "%02d", k);
    std::string base = ix_dir + "/test.00." + std::string(kk) + "mer";
    CHECK(kix.open(base + ".kix"));
    CHECK(kpx.open(base + ".kpx"));
    CHECK(ksx.open(base + ".ksx"));

    SearchConfig config;
    OidFilter no_filter;

    std::vector<SearchResult> local_results;
    for (const auto& q : queries) {
        auto sr = search_volume<uint16_t>(q.id, q.sequence, k,
                                          kix, kpx, ksx, no_filter, config);
        local_results.push_back(sr);
    }

    // --- Start ikafssnserver on UNIX socket ---
    std::string sock_path = g_test_dir + "/test_hc_server.sock";
    ::unlink(sock_path.c_str());

    ServerConfig sconfig;
    sconfig.ix_dir = ix_dir;
    sconfig.unix_socket_path = sock_path;
    sconfig.num_threads = 2;
    sconfig.log_level = Logger::kError;

    Server server;
    std::thread server_thread([&] {
        server.run(sconfig);
    });

    // Wait for server to start listening
    for (int i = 0; i < 20; i++) {
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
        struct stat st;
        if (stat(sock_path.c_str(), &st) == 0) break;
    }

    // --- Start ikafssnhttpd (Drogon) ---
    uint16_t http_port = 18923; // use a high port unlikely to conflict
    auto backend = std::make_shared<BackendClient>(
        BackendMode::kUnix, sock_path);
    HttpController controller(backend);
    controller.register_routes("");

    drogon::app()
        .addListener("127.0.0.1", http_port)
        .setThreadNum(1)
        .setLogLevel(trantor::Logger::kFatal);

    std::thread httpd_thread([] {
        drogon::app().run();
    });

    // Wait for Drogon to start
    std::this_thread::sleep_for(std::chrono::milliseconds(500));

    // --- Test http_search ---
    std::string http_url =
        "http://127.0.0.1:" + std::to_string(http_port);

    SearchRequest req;
    req.k = static_cast<uint8_t>(k);
    for (const auto& q : queries) {
        req.queries.push_back({q.id, q.sequence});
    }

    SearchResponse resp;
    std::string error_msg;
    bool ok2 = http_search(http_url, req, resp, error_msg);
    if (!ok2) {
        std::fprintf(stderr, "  http_search failed: %s\n", error_msg.c_str());
    }
    CHECK(ok2);

    // --- Verify results match local search ---
    CHECK(resp.status == 0);
    CHECK_EQ(resp.k, static_cast<uint8_t>(k));
    CHECK_EQ(resp.results.size(), local_results.size());

    for (size_t qi = 0;
         qi < resp.results.size() && qi < local_results.size(); qi++) {
        const auto& http_qr = resp.results[qi];
        const auto& local_sr = local_results[qi];
        CHECK(http_qr.query_id == local_sr.query_id);
        CHECK_EQ(http_qr.hits.size(), local_sr.hits.size());

        for (size_t hi = 0;
             hi < http_qr.hits.size() && hi < local_sr.hits.size(); hi++) {
            const auto& hh = http_qr.hits[hi];
            const auto& lh = local_sr.hits[hi];

            std::string local_acc(ksx.accession(lh.seq_id));
            CHECK(hh.accession == local_acc);
            CHECK_EQ(hh.q_start, lh.q_start);
            CHECK_EQ(hh.q_end, lh.q_end);
            CHECK_EQ(hh.s_start, lh.s_start);
            CHECK_EQ(hh.s_end, lh.s_end);
            CHECK_EQ(hh.score, static_cast<uint16_t>(lh.score));
            bool hh_reverse = (hh.strand == 1);
            CHECK(hh_reverse == lh.is_reverse);
        }
    }

    // --- Test seqidlist filtering through full stack ---
    std::string target_acc;
    if (ksx.num_sequences() > 0) {
        target_acc = std::string(ksx.accession(0));
    }
    if (!target_acc.empty()) {
        SearchRequest freq;
        freq.k = static_cast<uint8_t>(k);
        freq.seqidlist_mode = SeqidlistMode::kInclude;
        freq.seqids = {target_acc};
        for (const auto& q : queries) {
            freq.queries.push_back({q.id, q.sequence});
        }

        SearchResponse fresp;
        std::string ferr;
        CHECK(http_search(http_url, freq, fresp, ferr));
        CHECK(fresp.status == 0);

        for (const auto& qr : fresp.results) {
            for (const auto& hit : qr.hits) {
                CHECK(hit.accession == target_acc);
            }
        }
    }

    // --- Shutdown ---
    drogon::app().quit();
    httpd_thread.join();

    server.request_shutdown();
    server_thread.join();
    ::unlink(sock_path.c_str());
}

int main() {
    std::string source_dir = SOURCE_DIR;
    g_testdb_path = source_dir + "/test/testdata/testdb";

    struct stat st;
    if (stat((g_testdb_path + ".ndb").c_str(), &st) != 0 &&
        stat((g_testdb_path + ".nsq").c_str(), &st) != 0) {
        std::fprintf(stderr,
            "Test BLAST DB not found at %s. "
            "Run test/scripts/create_test_blastdb.sh first.\n",
            g_testdb_path.c_str());
        return 1;
    }

    g_test_dir = "/tmp/ikafssn_test_httpd_client_" +
                 std::to_string(::getpid());
    std::filesystem::create_directories(g_test_dir);

    test_http_client_search();

    std::filesystem::remove_all(g_test_dir);

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
