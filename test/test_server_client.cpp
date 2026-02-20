#include "test_util.hpp"

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

#include <cstdio>
#include <filesystem>
#include <string>
#include <thread>
#include <chrono>

#include <sys/stat.h>
#include <unistd.h>

using namespace ikafssn;

static std::string g_testdb_path;
static std::string g_test_dir;

// Build test index and return the index directory path
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
    bconfig.partitions = 1;
    bconfig.buffer_size = uint64_t(1) << 30;

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

    return ix_dir;
}

// Test: direct local search produces results, and server produces same results
static void test_server_client_search() {
    std::fprintf(stderr, "-- test_server_client_search\n");

    int k = 7;
    std::string ix_dir = build_test_index(k);
    CHECK(!ix_dir.empty());

    // Read query from test FASTA (use sequence from BLAST DB itself for guaranteed hits)
    std::string query_fasta = std::string(SOURCE_DIR) + "/test/testdata/queries.fasta";
    auto queries = read_fasta(query_fasta);
    CHECK(!queries.empty());

    // --- Direct local search for reference results ---
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

    // --- Set up server in a thread ---
    std::string sock_path = g_test_dir + "/test_server.sock";
    // Remove stale socket
    ::unlink(sock_path.c_str());

    // Load server index
    Server server;
    Logger logger(Logger::kError);
    CHECK(server.load_indexes(ix_dir, logger));
    CHECK(server.default_k() == k);

    // Start a listening socket
    int listen_fd = unix_listen(sock_path);
    CHECK(listen_fd >= 0);

    // Accept one connection in a background thread and process it
    std::thread server_thread([&] {
        int client_fd = accept_connection(listen_fd);
        if (client_fd >= 0) {
            handle_connection(client_fd, server.kmer_groups(),
                              server.default_k(), config, logger);
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

        for (size_t hi = 0; hi < server_qr.hits.size() && hi < local_sr.hits.size(); hi++) {
            const auto& sh = server_qr.hits[hi];
            const auto& lh = local_sr.hits[hi];

            std::string local_acc(ksx.accession(lh.seq_id));
            CHECK(sh.accession == local_acc);
            CHECK_EQ(sh.q_start, lh.q_start);
            CHECK_EQ(sh.q_end, lh.q_end);
            CHECK_EQ(sh.s_start, lh.s_start);
            CHECK_EQ(sh.s_end, lh.s_end);
            CHECK_EQ(sh.score, static_cast<uint16_t>(lh.score));
            bool sh_reverse = (sh.strand == 1);
            CHECK(sh_reverse == lh.is_reverse);
        }
    }
}

// Test: health check
static void test_health_check() {
    std::fprintf(stderr, "-- test_health_check\n");

    std::string sock_path = g_test_dir + "/test_health.sock";
    ::unlink(sock_path.c_str());

    int k = 7;
    std::string ix_dir = g_test_dir + "/sc_index";

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
    std::string ix_dir = g_test_dir + "/sc_index";
    std::string sock_path = g_test_dir + "/test_seqidlist.sock";
    ::unlink(sock_path.c_str());

    // Read queries
    std::string query_fasta = std::string(SOURCE_DIR) + "/test/testdata/queries.fasta";
    auto queries = read_fasta(query_fasta);
    CHECK(!queries.empty());

    // Get an accession from the ksx to use as seqidlist filter
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
    std::string source_dir = SOURCE_DIR;
    g_testdb_path = source_dir + "/test/testdata/testdb";

    // Check if test BLAST DB exists
    struct stat st;
    if (stat((g_testdb_path + ".ndb").c_str(), &st) != 0 &&
        stat((g_testdb_path + ".nsq").c_str(), &st) != 0) {
        std::fprintf(stderr, "Test BLAST DB not found at %s. "
                     "Run test/scripts/create_test_blastdb.sh first.\n",
                     g_testdb_path.c_str());
        return 1;
    }

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
