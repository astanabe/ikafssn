#include "ikafssnclient/socket_client.hpp"
#ifdef IKAFSSN_ENABLE_HTTP
#include "ikafssnclient/http_client.hpp"
#endif
#include "io/fasta_reader.hpp"
#include "io/seqidlist_reader.hpp"
#include "io/result_writer.hpp"
#include "protocol/messages.hpp"
#include "util/cli_parser.hpp"
#include "util/socket_utils.hpp"
#include "util/logger.hpp"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <thread>
#include <unordered_map>

using namespace ikafssn;

static void print_usage(const char* prog) {
    std::fprintf(stderr,
        "Usage: %s [options]\n"
        "\n"
        "Connection (one required):\n"
        "  -socket <path>           UNIX domain socket path\n"
        "  -tcp <host>:<port>       TCP server address\n"
#ifdef IKAFSSN_ENABLE_HTTP
        "  -http <url>              ikafssnhttpd URL (e.g., http://example.com:8080)\n"
#endif
        "\n"
        "Required:\n"
        "  -query <path>            Query FASTA file (- for stdin)\n"
        "\n"
        "Options:\n"
        "  -o <path>                Output file (default: stdout)\n"
        "  -k <int>                 K-mer size (default: server default)\n"
        "  -mode <1|2>              1=Stage1 only, 2=Stage1+Stage2 (default: server default)\n"
        "  -stage1_score <1|2>      1=coverscore, 2=matchscore (default: server default)\n"
        "  -sort_score <1|2>        1=stage1 score, 2=chainscore (default: server default)\n"
        "  -min_score <int>         Minimum score (default: server default)\n"
        "  -max_gap <int>           Chaining gap tolerance (default: server default)\n"
        "  -max_freq <int>          High-freq k-mer skip threshold (default: server default)\n"
        "  -min_diag_hits <int>     Diagonal filter min hits (default: server default)\n"
        "  -stage1_topn <int>       Stage 1 candidate limit (default: server default)\n"
        "  -min_stage1_score <num>  Stage 1 minimum score; integer or 0<P<1 fraction (default: server default)\n"
        "  -num_results <int>       Max results per query (default: server default)\n"
        "  -seqidlist <path>        Include only listed accessions\n"
        "  -negative_seqidlist <path>  Exclude listed accessions\n"
        "  -strand <-1|1|2>         Strand: 1=plus, -1=minus, 2=both (default: server default)\n"
        "  -accept_qdegen <0|1>     Accept queries with degenerate bases (default: 0)\n"
        "  -outfmt <tab|json>       Output format (default: tab)\n"
        "  -v, --verbose            Verbose logging\n"
#ifdef IKAFSSN_ENABLE_HTTP
        "\n"
        "HTTP Authentication:\n"
        "  --user <user:password>   Credentials (curl-style)\n"
        "  --http-user <USER>       Username (wget-style)\n"
        "  --http-password <PASS>   Password (wget-style, used with --http-user)\n"
        "  --netrc-file <path>      .netrc file for credentials\n"
#endif
        ,
        prog);
}

// Execute a search request via socket or HTTP, returning true on success.
static bool execute_search(
    const CliParser& cli,
    bool has_http,
    const SearchRequest& req,
    SearchResponse& resp,
    const Logger& logger
#ifdef IKAFSSN_ENABLE_HTTP
    , const HttpAuthConfig& auth
#endif
    ) {

#ifdef IKAFSSN_ENABLE_HTTP
    if (has_http) {
        std::string http_url = cli.get_string("-http");
        logger.debug("Connecting via HTTP to %s", http_url.c_str());

        std::string error_msg;
        if (!http_search(http_url, req, resp, error_msg, auth)) {
            std::fprintf(stderr, "Error: %s\n", error_msg.c_str());
            return false;
        }
        return true;
    }
#else
    (void)has_http;
#endif

    // Socket mode: create a new connection per attempt
    int fd = -1;
    if (cli.has("-socket")) {
        std::string sock_path = cli.get_string("-socket");
        fd = unix_connect(sock_path);
        if (fd < 0) {
            std::fprintf(stderr, "Error: cannot connect to UNIX socket %s\n", sock_path.c_str());
            return false;
        }
        logger.debug("Connected to UNIX socket %s", sock_path.c_str());
    } else {
        std::string tcp_addr = cli.get_string("-tcp");
        fd = tcp_connect(tcp_addr);
        if (fd < 0) {
            std::fprintf(stderr, "Error: cannot connect to TCP %s\n", tcp_addr.c_str());
            return false;
        }
        logger.debug("Connected to TCP %s", tcp_addr.c_str());
    }

    if (!socket_search(fd, req, resp)) {
        std::fprintf(stderr, "Error: search request failed\n");
        close_fd(fd);
        return false;
    }
    close_fd(fd);
    return true;
}

int main(int argc, char* argv[]) {
    CliParser cli(argc, argv);

    if (cli.has("-h") || cli.has("--help")) {
        print_usage(argv[0]);
        return 0;
    }

    bool has_http = false;
#ifdef IKAFSSN_ENABLE_HTTP
    has_http = cli.has("-http");
#endif

    if (!cli.has("-socket") && !cli.has("-tcp") && !has_http) {
        std::fprintf(stderr, "Error: one of -socket, -tcp"
#ifdef IKAFSSN_ENABLE_HTTP
            ", or -http"
#endif
            " is required\n");
        print_usage(argv[0]);
        return 1;
    }

    if (!cli.has("-query")) {
        std::fprintf(stderr, "Error: -query is required\n");
        print_usage(argv[0]);
        return 1;
    }

    // Check mutually exclusive options
    if (cli.has("-seqidlist") && cli.has("-negative_seqidlist")) {
        std::fprintf(stderr, "Error: -seqidlist and -negative_seqidlist are mutually exclusive\n");
        return 1;
    }

    bool verbose = cli.has("-v") || cli.has("--verbose");
    Logger logger(verbose ? Logger::kDebug : Logger::kInfo);

#ifdef IKAFSSN_ENABLE_HTTP
    // HTTP authentication resolution
    HttpAuthConfig auth;
    if (cli.has("--user") && cli.has("--http-user")) {
        std::fprintf(stderr, "Error: --user and --http-user are mutually exclusive\n");
        return 1;
    }
    if (cli.has("--user")) {
        auth.userpwd = cli.get_string("--user");
    } else if (cli.has("--http-user")) {
        std::string user = cli.get_string("--http-user");
        std::string pass = cli.get_string("--http-password", "");
        auth.userpwd = user + ":" + pass;
    }
    if (cli.has("--netrc-file")) {
        auth.netrc_file = cli.get_string("--netrc-file");
    }
#endif

    std::string query_path = cli.get_string("-query");
    std::string output_path = cli.get_string("-o");

    // Output format
    OutputFormat outfmt = OutputFormat::kTab;
    std::string outfmt_str = cli.get_string("-outfmt", "tab");
    if (outfmt_str == "json") {
        outfmt = OutputFormat::kJson;
    } else if (outfmt_str != "tab") {
        std::fprintf(stderr, "Error: unknown output format '%s'\n", outfmt_str.c_str());
        return 1;
    }

    // Read query FASTA
    auto queries = read_fasta(query_path);
    if (queries.empty()) {
        std::fprintf(stderr, "Error: no query sequences found\n");
        return 1;
    }
    logger.info("Read %zu query sequence(s)", queries.size());

    // Build base search request parameters (shared across retries)
    SearchRequest base_req;
    base_req.k = static_cast<uint8_t>(cli.get_int("-k", 0));
    base_req.min_score = static_cast<uint16_t>(cli.get_int("-min_score", 0));
    base_req.max_gap = static_cast<uint16_t>(cli.get_int("-max_gap", 0));
    base_req.max_freq = static_cast<uint32_t>(cli.get_int("-max_freq", 0));
    base_req.min_diag_hits = static_cast<uint8_t>(cli.get_int("-min_diag_hits", 0));
    base_req.stage1_topn = static_cast<uint16_t>(cli.get_int("-stage1_topn", 0));
    {
        double min_s1 = cli.get_double("-min_stage1_score", 0.0);
        if (min_s1 > 0 && min_s1 < 1.0) {
            base_req.min_stage1_score_frac_x10000 =
                static_cast<uint16_t>(min_s1 * 10000.0);
        } else {
            base_req.min_stage1_score = static_cast<uint16_t>(min_s1);
        }
    }
    base_req.num_results = static_cast<uint16_t>(cli.get_int("-num_results", 0));
    base_req.mode = static_cast<uint8_t>(cli.get_int("-mode", 0));
    base_req.stage1_score_type = static_cast<uint8_t>(cli.get_int("-stage1_score", 0));
    base_req.sort_score = static_cast<uint8_t>(cli.get_int("-sort_score", 0));
    base_req.accept_qdegen = static_cast<uint8_t>(cli.get_int("-accept_qdegen", 0));
    base_req.strand = static_cast<int8_t>(cli.get_int("-strand", 0));

    // Seqidlist
    if (cli.has("-seqidlist")) {
        base_req.seqidlist_mode = SeqidlistMode::kInclude;
        base_req.seqids = read_seqidlist(cli.get_string("-seqidlist"));
        logger.info("Loaded %zu accessions (include mode)", base_req.seqids.size());
    } else if (cli.has("-negative_seqidlist")) {
        base_req.seqidlist_mode = SeqidlistMode::kExclude;
        base_req.seqids = read_seqidlist(cli.get_string("-negative_seqidlist"));
        logger.info("Loaded %zu accessions (exclude mode)", base_req.seqids.size());
    }

    // Build query_id -> sequence lookup for retry
    std::unordered_map<std::string, std::string> query_map;
    for (const auto& q : queries) {
        query_map[q.id] = q.sequence;
    }

    // Build initial request with all queries
    SearchRequest req = base_req;
    for (const auto& q : queries) {
        req.queries.push_back({q.id, q.sequence});
    }

    // Accumulate results across retry attempts
    std::vector<OutputHit> all_hits;
    bool has_skipped = false;
    uint8_t resp_mode = 0;
    uint8_t resp_stage1_score_type = 0;

    // Retry schedule: 30s, 60s, 120s, 120s, ...
    static constexpr int retry_delays[] = {30, 60, 120};
    static constexpr int num_retry_delays =
        static_cast<int>(sizeof(retry_delays) / sizeof(retry_delays[0]));

    for (int attempt = 0; ; attempt++) {
        SearchResponse resp;
        if (!execute_search(cli, has_http, req, resp, logger
#ifdef IKAFSSN_ENABLE_HTTP
                , auth
#endif
                )) {
            return 1;
        }

        if (resp.status != 0) {
            std::fprintf(stderr, "Error: server returned status %d\n", resp.status);
            return 1;
        }

        logger.info("Received response: k=%d, %zu query result(s), %zu rejected",
                     resp.k, resp.results.size(), resp.rejected_query_ids.size());

        // Save mode/score_type from first successful response
        if (attempt == 0) {
            resp_mode = resp.mode;
            resp_stage1_score_type = resp.stage1_score_type;
        }

        // Collect results from accepted queries
        for (const auto& qr : resp.results) {
            if (qr.skipped != 0) {
                has_skipped = true;
                std::fprintf(stderr, "Warning: query '%s' was skipped (degenerate bases)\n",
                             qr.query_id.c_str());
                continue;
            }
            for (const auto& hit : qr.hits) {
                OutputHit oh;
                oh.query_id = qr.query_id;
                oh.accession = hit.accession;
                oh.strand = (hit.strand == 0) ? '+' : '-';
                oh.q_start = hit.q_start;
                oh.q_end = hit.q_end;
                oh.s_start = hit.s_start;
                oh.s_end = hit.s_end;
                oh.score = hit.score;
                oh.stage1_score = hit.stage1_score;
                oh.volume = hit.volume;
                all_hits.push_back(oh);
            }
        }

        // If no rejected queries, we're done
        if (resp.rejected_query_ids.empty()) break;

        // Sleep before retry
        int delay_idx = std::min(attempt, num_retry_delays - 1);
        int delay = retry_delays[delay_idx];
        logger.info("%zu queries rejected, retrying in %d seconds...",
                    resp.rejected_query_ids.size(), delay);
        std::this_thread::sleep_for(std::chrono::seconds(delay));

        // Rebuild request with only rejected query IDs
        req = base_req;
        for (const auto& qid : resp.rejected_query_ids) {
            auto it = query_map.find(qid);
            if (it != query_map.end()) {
                req.queries.push_back({qid, it->second});
            }
        }

        if (req.queries.empty()) break; // all resolved
    }

    // Write output using mode/stage1_score_type from response
    if (output_path.empty()) {
        write_results(std::cout, all_hits, outfmt,
                      resp_mode, resp_stage1_score_type);
    } else {
        std::ofstream out(output_path);
        if (!out.is_open()) {
            std::fprintf(stderr, "Error: cannot open output file %s\n", output_path.c_str());
            return 1;
        }
        write_results(out, all_hits, outfmt,
                      resp_mode, resp_stage1_score_type);
    }

    logger.info("Done. %zu hit(s) reported.", all_hits.size());
    return has_skipped ? 2 : 0;
}
