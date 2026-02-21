#include "test_util.hpp"
#include "ssu_test_fixture.hpp"

#include "io/blastdb_reader.hpp"
#include "io/fasta_reader.hpp"
#include "io/result_writer.hpp"
#include "index/index_builder.hpp"
#include "index/kix_reader.hpp"
#include "index/kpx_reader.hpp"
#include "index/ksx_reader.hpp"
#include "search/oid_filter.hpp"
#include "search/volume_searcher.hpp"
#include "ikafssnserver/server.hpp"
#include "ikafssnserver/connection_handler.hpp"
#include "ikafssnserver/request_processor.hpp"
#include "ikafssnclient/socket_client.hpp"
#include "protocol/frame.hpp"
#include "protocol/messages.hpp"
#include "protocol/serializer.hpp"
#include "util/socket_utils.hpp"
#include "util/logger.hpp"
#include "core/config.hpp"

#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <string>
#include <thread>
#include <chrono>

#include <tbb/task_arena.h>

#include <sys/stat.h>
#include <unistd.h>

using namespace ikafssn;
using namespace ssu_fixture;

static std::string g_testdb_path;
static std::string g_test_dir;

// Build test index and return the index prefix path
static std::string build_test_index(int k) {
    std::string ix_dir = g_test_dir + "/sc_index";
    std::filesystem::create_directories(ix_dir);

    BlastDbReader db;
    if (!db.open(g_testdb_path)) {
        std::fprintf(stderr, "Cannot open test BLAST DB at %s\n", g_testdb_path.c_str());
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

    if (!ok) {
        std::fprintf(stderr, "Failed to build index for k=%d\n", k);
        return {};
    }

    return ix_dir + "/test";
}

// Test: direct local search produces results, and server produces same results
static void test_server_client_search() {
    std::fprintf(stderr, "-- test_server_client_search\n");

    int k = 7;
    std::string ix_prefix = build_test_index(k);
    CHECK(!ix_prefix.empty());

    // Read query from derived test data
    std::string query_fasta = queries_path();
    auto queries = read_fasta(query_fasta);
    CHECK(!queries.empty());

    // --- Direct local search for reference results ---
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
    OidFilter no_filter;

    std::vector<SearchResult> local_results;
    for (const auto& q : queries) {
        auto sr = search_volume<uint16_t>(q.id, q.sequence, k,
                                          kix, kpx, ksx, no_filter, config);
        local_results.push_back(sr);
    }

    // --- Set up server in a thread ---
    std::string sock_path = g_test_dir + "/test_server.sock";
    // Remove stale socket
    ::unlink(sock_path.c_str());

    // Load server index
    Server server;
    Logger logger(Logger::kError);
    CHECK(server.load_indexes(ix_prefix, logger));
    CHECK(server.default_k() == k);

    // Start a listening socket
    int listen_fd = unix_listen(sock_path);
    CHECK(listen_fd >= 0);

    // Accept one connection in a background thread and process it
    tbb::task_arena arena(1);
    std::thread server_thread([&] {
        int client_fd = accept_connection(listen_fd);
        if (client_fd >= 0) {
            handle_connection(client_fd, server.kmer_groups(),
                              server.default_k(), config, server, arena, logger);
        }
    });

    // Give server thread time to reach accept()
    std::this_thread::sleep_for(std::chrono::milliseconds(50));

    // --- Client connects and searches ---
    int fd = unix_connect(sock_path);
    CHECK(fd >= 0);

    SearchRequest req;
    req.k = static_cast<uint8_t>(k);
    for (const auto& q : queries) {
        req.queries.push_back({q.id, q.sequence});
    }

    SearchResponse resp;
    CHECK(socket_search(fd, req, resp));
    close_fd(fd);

    server_thread.join();
    close_fd(listen_fd);
    ::unlink(sock_path.c_str());

    // --- Verify results match ---
    CHECK(resp.status == 0);
    CHECK_EQ(resp.k, static_cast<uint8_t>(k));
    CHECK_EQ(resp.results.size(), local_results.size());

    for (size_t qi = 0; qi < resp.results.size() && qi < local_results.size(); qi++) {
        const auto& server_qr = resp.results[qi];
        const auto& local_sr = local_results[qi];
        CHECK(server_qr.query_id == local_sr.query_id);
        CHECK_EQ(server_qr.hits.size(), local_sr.hits.size());

        // Build comparable lists sorted by (accession, s_start, q_start) to avoid ordering issues
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

        std::vector<HitKey> server_sorted, local_sorted;
        for (const auto& sh : server_qr.hits) {
            server_sorted.push_back({sh.accession, sh.q_start, sh.q_end,
                                     sh.s_start, sh.s_end, sh.score, sh.strand == 1});
        }
        for (const auto& lh : local_sr.hits) {
            local_sorted.push_back({std::string(ksx.accession(lh.seq_id)),
                                    lh.q_start, lh.q_end, lh.s_start, lh.s_end,
                                    static_cast<uint16_t>(lh.score), lh.is_reverse});
        }
        std::sort(server_sorted.begin(), server_sorted.end());
        std::sort(local_sorted.begin(), local_sorted.end());

        for (size_t hi = 0; hi < server_sorted.size() && hi < local_sorted.size(); hi++) {
            CHECK(server_sorted[hi].accession == local_sorted[hi].accession);
            CHECK_EQ(server_sorted[hi].q_start, local_sorted[hi].q_start);
            CHECK_EQ(server_sorted[hi].q_end, local_sorted[hi].q_end);
            CHECK_EQ(server_sorted[hi].s_start, local_sorted[hi].s_start);
            CHECK_EQ(server_sorted[hi].s_end, local_sorted[hi].s_end);
            CHECK_EQ(server_sorted[hi].score, local_sorted[hi].score);
            CHECK(server_sorted[hi].is_reverse == local_sorted[hi].is_reverse);
        }
    }
}

// Test: health check
static void test_health_check() {
    std::fprintf(stderr, "-- test_health_check\n");

    std::string sock_path = g_test_dir + "/test_health.sock";
    ::unlink(sock_path.c_str());

    int k = 7;
    std::string ix_prefix = g_test_dir + "/sc_index/test";

    Server server;
    Logger logger(Logger::kError);
    CHECK(server.load_indexes(ix_prefix, logger));

    SearchConfig config;
    int listen_fd = unix_listen(sock_path);
    CHECK(listen_fd >= 0);

    tbb::task_arena arena(1);
    std::thread server_thread([&] {
        int client_fd = accept_connection(listen_fd);
        if (client_fd >= 0) {
            handle_connection(client_fd, server.kmer_groups(),
                              server.default_k(), config, server, arena, logger);
        }
    });

    std::this_thread::sleep_for(std::chrono::milliseconds(50));

    int fd = unix_connect(sock_path);
    CHECK(fd >= 0);

    HealthResponse hresp;
    CHECK(socket_health_check(fd, hresp));
    CHECK_EQ(hresp.status, 0);

    close_fd(fd);
    server_thread.join();
    close_fd(listen_fd);
    ::unlink(sock_path.c_str());
}

// Test: seqidlist filtering through server
static void test_seqidlist_filter_via_server() {
    std::fprintf(stderr, "-- test_seqidlist_filter_via_server\n");

    int k = 7;
    std::string ix_prefix = g_test_dir + "/sc_index/test";
    std::string sock_path = g_test_dir + "/test_seqidlist.sock";
    ::unlink(sock_path.c_str());

    // Read queries
    auto queries = read_fasta(queries_path());
    CHECK(!queries.empty());

    // Get an accession from the ksx to use as seqidlist filter
    KsxReader ksx;
    char kk[4];
    std::snprintf(kk, sizeof(kk), "%02d", k);
    std::string base = ix_prefix + ".00." + std::string(kk) + "mer";
    CHECK(ksx.open(base + ".ksx"));

    std::string target_acc;
    if (ksx.num_sequences() > 0) {
        target_acc = std::string(ksx.accession(0));
    }
    CHECK(!target_acc.empty());

    Server server;
    Logger logger(Logger::kError);
    CHECK(server.load_indexes(ix_prefix, logger));

    SearchConfig config;
    int listen_fd = unix_listen(sock_path);
    CHECK(listen_fd >= 0);

    tbb::task_arena arena(1);
    std::thread server_thread([&] {
        int client_fd = accept_connection(listen_fd);
        if (client_fd >= 0) {
            handle_connection(client_fd, server.kmer_groups(),
                              server.default_k(), config, server, arena, logger);
        }
    });

    std::this_thread::sleep_for(std::chrono::milliseconds(50));

    int fd = unix_connect(sock_path);
    CHECK(fd >= 0);

    // Build request with seqidlist (include mode)
    SearchRequest req;
    req.k = static_cast<uint8_t>(k);
    req.seqidlist_mode = SeqidlistMode::kInclude;
    req.seqids = {target_acc};
    for (const auto& q : queries) {
        req.queries.push_back({q.id, q.sequence});
    }

    SearchResponse resp;
    CHECK(socket_search(fd, req, resp));
    close_fd(fd);

    server_thread.join();
    close_fd(listen_fd);
    ::unlink(sock_path.c_str());

    CHECK(resp.status == 0);

    // All hits should have the target accession
    for (const auto& qr : resp.results) {
        for (const auto& hit : qr.hits) {
            CHECK(hit.accession == target_acc);
        }
    }
}

int main() {
    check_ssu_available();
    check_derived_data_ready();

    g_testdb_path = ssu_db_prefix();

    // Create temp directory
    g_test_dir = "/tmp/ikafssn_test_server_client_" + std::to_string(::getpid());
    std::filesystem::create_directories(g_test_dir);

    test_server_client_search();
    test_health_check();
    test_seqidlist_filter_via_server();

    // Cleanup
    std::filesystem::remove_all(g_test_dir);

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
