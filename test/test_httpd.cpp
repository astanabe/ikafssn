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
#include "ikafssnserver/connection_handler.hpp"
#include "ikafssnhttpd/backend_client.hpp"
#include "protocol/messages.hpp"
#include "util/socket_utils.hpp"
#include "util/logger.hpp"
#include "core/config.hpp"

#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <string>
#include <thread>
#include <chrono>

#include <sys/stat.h>
#include <unistd.h>

using namespace ikafssn;
using namespace ssu_fixture;

static std::string g_testdb_path;
static std::string g_test_dir;

static std::string build_test_index(int k) {
    std::string ix_dir = g_test_dir + "/httpd_index";
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

// Test: BackendClient search matches direct local search
static void test_backend_search() {
    std::fprintf(stderr, "-- test_backend_search\n");

    int k = 7;
    std::string ix_dir = build_test_index(k);
    CHECK(!ix_dir.empty());

    // Read query FASTA
    auto queries = read_fasta(queries_path());
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

    // --- Start server on UNIX socket ---
    std::string sock_path = g_test_dir + "/test_httpd_search.sock";
    ::unlink(sock_path.c_str());

    Server server;
    Logger logger(Logger::kError);
    CHECK(server.load_indexes(ix_dir, logger));

    int listen_fd = unix_listen(sock_path);
    CHECK(listen_fd >= 0);

    // Accept one connection and handle it
    std::thread server_thread([&] {
        int client_fd = accept_connection(listen_fd);
        if (client_fd >= 0) {
            handle_connection(client_fd, server.kmer_groups(),
                              server.default_k(), config, logger);
        }
    });

    std::this_thread::sleep_for(std::chrono::milliseconds(50));

    // --- Search via BackendClient ---
    BackendClient backend(BackendMode::kUnix, sock_path);

    SearchRequest req;
    req.k = static_cast<uint8_t>(k);
    for (const auto& q : queries) {
        req.queries.push_back({q.id, q.sequence});
    }

    SearchResponse resp;
    std::string error_msg;
    CHECK(backend.search(req, resp, error_msg));

    server_thread.join();
    close_fd(listen_fd);
    ::unlink(sock_path.c_str());

    // --- Verify results match ---
    CHECK(resp.status == 0);
    CHECK_EQ(resp.k, static_cast<uint8_t>(k));
    CHECK_EQ(resp.results.size(), local_results.size());

    for (size_t qi = 0;
         qi < resp.results.size() && qi < local_results.size(); qi++) {
        const auto& server_qr = resp.results[qi];
        const auto& local_sr = local_results[qi];
        CHECK(server_qr.query_id == local_sr.query_id);
        CHECK_EQ(server_qr.hits.size(), local_sr.hits.size());

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

        for (size_t hi = 0;
             hi < server_sorted.size() && hi < local_sorted.size(); hi++) {
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

// Test: BackendClient health check
static void test_backend_health() {
    std::fprintf(stderr, "-- test_backend_health\n");

    std::string sock_path = g_test_dir + "/test_httpd_health.sock";
    ::unlink(sock_path.c_str());

    std::string ix_dir = g_test_dir + "/httpd_index";

    Server server;
    Logger logger(Logger::kError);
    CHECK(server.load_indexes(ix_dir, logger));

    SearchConfig config;
    int listen_fd = unix_listen(sock_path);
    CHECK(listen_fd >= 0);

    std::thread server_thread([&] {
        int client_fd = accept_connection(listen_fd);
        if (client_fd >= 0) {
            handle_connection(client_fd, server.kmer_groups(),
                              server.default_k(), config, logger);
        }
    });

    std::this_thread::sleep_for(std::chrono::milliseconds(50));

    BackendClient backend(BackendMode::kUnix, sock_path);

    HealthResponse hresp;
    std::string error_msg;
    CHECK(backend.health_check(hresp, error_msg));
    CHECK_EQ(hresp.status, 0);

    server_thread.join();
    close_fd(listen_fd);
    ::unlink(sock_path.c_str());
}

// Test: BackendClient info request
static void test_backend_info() {
    std::fprintf(stderr, "-- test_backend_info\n");

    std::string sock_path = g_test_dir + "/test_httpd_info.sock";
    ::unlink(sock_path.c_str());

    std::string ix_dir = g_test_dir + "/httpd_index";

    Server server;
    Logger logger(Logger::kError);
    CHECK(server.load_indexes(ix_dir, logger));

    SearchConfig config;
    int listen_fd = unix_listen(sock_path);
    CHECK(listen_fd >= 0);

    std::thread server_thread([&] {
        int client_fd = accept_connection(listen_fd);
        if (client_fd >= 0) {
            handle_connection(client_fd, server.kmer_groups(),
                              server.default_k(), config, logger);
        }
    });

    std::this_thread::sleep_for(std::chrono::milliseconds(50));

    BackendClient backend(BackendMode::kUnix, sock_path);

    InfoResponse iresp;
    std::string error_msg;
    CHECK(backend.info(iresp, error_msg));
    CHECK_EQ(iresp.status, 0);
    CHECK_EQ(iresp.default_k, static_cast<uint8_t>(server.default_k()));

    // Should have at least one kmer group
    CHECK(!iresp.groups.empty());

    // First group should have k=7 (we built with k=7)
    CHECK_EQ(iresp.groups[0].k, 7);
    CHECK_EQ(iresp.groups[0].kmer_type, 0); // uint16 for k<9
    CHECK(!iresp.groups[0].volumes.empty());

    // Volume 0 should have 5 sequences (our test db)
    CHECK_EQ(iresp.groups[0].volumes[0].volume_index, 0);
    CHECK(iresp.groups[0].volumes[0].num_sequences > 0);
    CHECK(iresp.groups[0].volumes[0].total_postings > 0);

    server_thread.join();
    close_fd(listen_fd);
    ::unlink(sock_path.c_str());
}

// Test: BackendClient error handling (connection failure)
static void test_backend_connection_failure() {
    std::fprintf(stderr, "-- test_backend_connection_failure\n");

    // Connect to a non-existent socket
    BackendClient backend(BackendMode::kUnix,
                          g_test_dir + "/nonexistent.sock");

    SearchRequest req;
    req.queries.push_back({"q1", "ACGTACGT"});

    SearchResponse resp;
    std::string error_msg;
    CHECK(!backend.search(req, resp, error_msg));
    CHECK(!error_msg.empty());

    HealthResponse hresp;
    std::string health_error;
    CHECK(!backend.health_check(hresp, health_error));
    CHECK(!health_error.empty());
}

// Test: BackendClient seqidlist filtering
static void test_backend_seqidlist_filter() {
    std::fprintf(stderr, "-- test_backend_seqidlist_filter\n");

    int k = 7;
    std::string ix_dir = g_test_dir + "/httpd_index";
    std::string sock_path = g_test_dir + "/test_httpd_seqidlist.sock";
    ::unlink(sock_path.c_str());

    // Read queries
    auto queries = read_fasta(queries_path());
    CHECK(!queries.empty());

    // Get an accession from ksx
    KsxReader ksx;
    char kk[4];
    std::snprintf(kk, sizeof(kk), "%02d", k);
    std::string base = ix_dir + "/test.00." + std::string(kk) + "mer";
    CHECK(ksx.open(base + ".ksx"));

    std::string target_acc;
    if (ksx.num_sequences() > 0) {
        target_acc = std::string(ksx.accession(0));
    }
    CHECK(!target_acc.empty());

    Server server;
    Logger logger(Logger::kError);
    CHECK(server.load_indexes(ix_dir, logger));

    SearchConfig config;
    int listen_fd = unix_listen(sock_path);
    CHECK(listen_fd >= 0);

    std::thread server_thread([&] {
        int client_fd = accept_connection(listen_fd);
        if (client_fd >= 0) {
            handle_connection(client_fd, server.kmer_groups(),
                              server.default_k(), config, logger);
        }
    });

    std::this_thread::sleep_for(std::chrono::milliseconds(50));

    BackendClient backend(BackendMode::kUnix, sock_path);

    SearchRequest req;
    req.k = static_cast<uint8_t>(k);
    req.seqidlist_mode = SeqidlistMode::kInclude;
    req.seqids = {target_acc};
    for (const auto& q : queries) {
        req.queries.push_back({q.id, q.sequence});
    }

    SearchResponse resp;
    std::string error_msg;
    CHECK(backend.search(req, resp, error_msg));

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

    g_test_dir = "/tmp/ikafssn_test_httpd_" + std::to_string(::getpid());
    std::filesystem::create_directories(g_test_dir);

    test_backend_search();
    test_backend_health();
    test_backend_info();
    test_backend_connection_failure();
    test_backend_seqidlist_filter();

    // Cleanup
    std::filesystem::remove_all(g_test_dir);

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
