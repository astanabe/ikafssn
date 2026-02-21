#include "test_util.hpp"
#include "ssu_test_fixture.hpp"

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

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <string>
#include <thread>
#include <chrono>
#include <atomic>

#include <sys/stat.h>
#include <unistd.h>

using namespace ikafssn;
using namespace ssu_fixture;

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

    return ix_dir + "/test";
}

// Full stack test: ikafssnclient (http_search) -> ikafssnhttpd -> ikafssnserver
static void test_http_client_search() {
    std::fprintf(stderr, "-- test_http_client_search\n");

    int k = 7;
    std::string ix_prefix = build_test_index(k);
    CHECK(!ix_prefix.empty());

    // Read query FASTA
    auto queries = read_fasta(queries_path());
    CHECK(!queries.empty());

    // --- Direct local search for reference ---
    KixReader kix;
    KpxReader kpx;
    KsxReader ksx;
    char kk[4];
    std::snprintf(kk, sizeof(kk), "%02d", k);
    std::string base = ix_prefix + ".00." + std::string(kk) + "mer";
    CHECK(kix.open(base + ".kix"));
    CHECK(kpx.open(base + ".kpx"));
    CHECK(ksx.open(base + ".ksx"));

    // Build search config matching server defaults.
    // ServerConfig.max_freq_raw defaults to 0.5 (fraction), which the server
    // resolves to ceil(0.5 * total_nseq).  Use the same resolved value for
    // the local reference search so results are comparable.
    SearchConfig config;
    {
        uint32_t total_nseq = ksx.num_sequences();
        auto resolved = static_cast<uint32_t>(
            std::ceil(0.5 * total_nseq));
        if (resolved == 0) resolved = 1;
        config.stage1.max_freq = resolved;
    }
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
    sconfig.ix_prefix = ix_prefix;
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

        // Sort both result sets before comparing to avoid ordering issues
        struct HitKey {
            std::string accession;
            uint32_t q_start, q_end, s_start, s_end;
            uint16_t score;
            bool is_reverse;
            bool operator<(const HitKey& o) const {
                if (accession != o.accession) return accession < o.accession;
                if (s_start != o.s_start) return s_start < o.s_start;
                return q_start < o.q_start;
            }
        };

        std::vector<HitKey> http_sorted, local_sorted;
        for (const auto& hh : http_qr.hits) {
            http_sorted.push_back({hh.accession, hh.q_start, hh.q_end,
                                   hh.s_start, hh.s_end, hh.score, hh.strand == 1});
        }
        for (const auto& lh : local_sr.hits) {
            local_sorted.push_back({std::string(ksx.accession(lh.seq_id)),
                                    lh.q_start, lh.q_end, lh.s_start, lh.s_end,
                                    static_cast<uint16_t>(lh.score), lh.is_reverse});
        }
        std::sort(http_sorted.begin(), http_sorted.end());
        std::sort(local_sorted.begin(), local_sorted.end());

        for (size_t hi = 0;
             hi < http_sorted.size() && hi < local_sorted.size(); hi++) {
            CHECK(http_sorted[hi].accession == local_sorted[hi].accession);
            CHECK_EQ(http_sorted[hi].q_start, local_sorted[hi].q_start);
            CHECK_EQ(http_sorted[hi].q_end, local_sorted[hi].q_end);
            CHECK_EQ(http_sorted[hi].s_start, local_sorted[hi].s_start);
            CHECK_EQ(http_sorted[hi].s_end, local_sorted[hi].s_end);
            CHECK_EQ(http_sorted[hi].score, local_sorted[hi].score);
            CHECK(http_sorted[hi].is_reverse == local_sorted[hi].is_reverse);
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
    check_ssu_available();
    check_derived_data_ready();

    g_testdb_path = ssu_db_prefix();

    g_test_dir = "/tmp/ikafssn_test_httpd_client_" +
                 std::to_string(::getpid());
    std::filesystem::create_directories(g_test_dir);

    test_http_client_search();

    std::filesystem::remove_all(g_test_dir);

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
