#include "ikafssnclient/socket_client.hpp"
#ifdef IKAFSSN_ENABLE_HTTP
#include "ikafssnclient/http_client.hpp"
#endif
#include "core/version.hpp"
#include "protocol/info_format.hpp"
#include "util/common_init.hpp"
#include "io/fasta_reader.hpp"
#include "io/seqidlist_reader.hpp"
#include "io/result_writer.hpp"
#include "io/sam_writer.hpp"
#include "protocol/messages.hpp"
#include "util/cli_parser.hpp"
#include "util/socket_utils.hpp"
#include "util/logger.hpp"

#include <chrono>
#include <climits>
#include <cstdio>
#include <cstdlib>
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
        "  -mode <1|2|3>            1=Stage1, 2=Stage1+2, 3=Stage1+2+3 (default: server default)\n"
        "  -stage1_score <1|2>      1=coverscore, 2=matchscore (default: server default)\n"
        "  -stage2_min_score <int>  Minimum chain score (default: server default)\n"
        "  -stage2_max_gap <int>    Chaining gap tolerance (default: server default)\n"
        "  -stage2_max_lookback <int>  Chaining DP lookback window (default: server default)\n"
        "  -stage1_max_freq <num>   High-freq k-mer skip threshold (default: server default)\n"
        "                           0 < x < 1: fraction of total NSEQ across all volumes\n"
        "                           >= 1: absolute count threshold\n"
        "  -stage2_min_diag_hits <int>  Diagonal filter min hits (default: server default)\n"
        "  -stage1_topn <int>       Stage 1 candidate limit (default: server default)\n"
        "  -stage1_min_score <num>  Stage 1 minimum score; integer or 0<P<1 fraction (default: server default)\n"
        "  -num_results <int>       Max results per query (default: server default)\n"
        "  -seqidlist <path>        Include only listed accessions\n"
        "  -negative_seqidlist <path>  Exclude listed accessions\n"
        "  -strand <-1|1|2>         Strand: 1=plus, -1=minus, 2=both (default: server default)\n"
        "  -db <name>               Target database name on server (required for multi-DB servers)\n"
        "  -accept_qdegen <0|1>     Accept queries with degenerate bases (default: 1)\n"
        "  -context <value>         Context extension (int=bases, decimal=ratio, default: 0)\n"
        "  -stage3_traceback <0|1>  Enable traceback in mode 3 (default: 0)\n"
        "  -stage3_gapopen <int>    Gap open penalty (default: server default)\n"
        "  -stage3_gapext <int>     Gap extension penalty (default: server default)\n"
        "  -stage3_min_pident <num> Min percent identity filter (default: server default)\n"
        "  -stage3_min_nident <int> Min identical bases filter (default: server default)\n"
        "  -outfmt <tab|json|sam|bam>  Output format (default: tab)\n"
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

// Execute an info request via socket or HTTP, returning true on success.
static bool execute_info(
    const CliParser& cli,
    bool has_http,
    InfoResponse& resp,
    const Logger& logger
#ifdef IKAFSSN_ENABLE_HTTP
    , const HttpAuthConfig& auth
#endif
    ) {

#ifdef IKAFSSN_ENABLE_HTTP
    if (has_http) {
        std::string http_url = cli.get_string("-http");
        logger.debug("Fetching server info via HTTP from %s", http_url.c_str());

        std::string error_msg;
        if (!http_info(http_url, resp, error_msg, auth)) {
            std::fprintf(stderr, "Error: failed to fetch server info: %s\n",
                         error_msg.c_str());
            return false;
        }
        return true;
    }
#else
    (void)has_http;
#endif

    int fd = -1;
    if (cli.has("-socket")) {
        std::string sock_path = cli.get_string("-socket");
        fd = unix_connect(sock_path);
        if (fd < 0) {
            std::fprintf(stderr, "Error: cannot connect to UNIX socket %s\n",
                         sock_path.c_str());
            return false;
        }
    } else {
        std::string tcp_addr = cli.get_string("-tcp");
        fd = tcp_connect(tcp_addr);
        if (fd < 0) {
            std::fprintf(stderr, "Error: cannot connect to TCP %s\n",
                         tcp_addr.c_str());
            return false;
        }
    }

    if (!socket_info(fd, resp)) {
        std::fprintf(stderr, "Error: info request failed\n");
        close_fd(fd);
        return false;
    }
    close_fd(fd);
    return true;
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

    if (check_version(cli, "ikafssnclient")) return 0;

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

    Logger logger = make_logger(cli);

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
    OutputFormat outfmt;
    {
        std::string err;
        if (!parse_output_format(cli.get_string("-outfmt", "tab"), outfmt, err)) {
            std::fprintf(stderr, "%s\n", err.c_str());
            return 1;
        }
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
    if (cli.has("-stage2_min_score")) {
        base_req.stage2_min_score = static_cast<uint16_t>(cli.get_int("-stage2_min_score", 0));
        base_req.has_stage2_min_score = 1;
    }
    base_req.stage2_max_gap = static_cast<uint16_t>(cli.get_int("-stage2_max_gap", 0));
    base_req.stage2_max_lookback = static_cast<uint16_t>(cli.get_int("-stage2_max_lookback", 0));
    {
        double max_freq_val = cli.get_double("-stage1_max_freq", 0.0);
        if (max_freq_val > 0 && max_freq_val < 1.0) {
            base_req.stage1_max_freq_frac_x10000 =
                static_cast<uint16_t>(max_freq_val * 10000.0);
        } else {
            base_req.stage1_max_freq = static_cast<uint32_t>(max_freq_val);
        }
    }
    base_req.stage2_min_diag_hits = static_cast<uint8_t>(cli.get_int("-stage2_min_diag_hits", 0));
    base_req.stage1_topn = static_cast<uint16_t>(cli.get_int("-stage1_topn", 0));
    {
        double min_s1 = cli.get_double("-stage1_min_score", 0.0);
        if (min_s1 > 0 && min_s1 < 1.0) {
            base_req.stage1_min_score_frac_x10000 =
                static_cast<uint16_t>(min_s1 * 10000.0);
        } else {
            base_req.stage1_min_score = static_cast<uint16_t>(min_s1);
        }
    }
    base_req.num_results = static_cast<uint16_t>(cli.get_int("-num_results", 0));
    base_req.mode = static_cast<uint8_t>(cli.get_int("-mode", 0));
    base_req.stage1_score = static_cast<uint8_t>(cli.get_int("-stage1_score", 0));
    base_req.accept_qdegen = static_cast<uint8_t>(cli.get_int("-accept_qdegen", 1));
    base_req.strand = static_cast<int8_t>(cli.get_int("-strand", 0));
    base_req.db = cli.get_string("-db", "");

    // Stage 3 parameters
    base_req.stage3_traceback = static_cast<uint8_t>(cli.get_int("-stage3_traceback", 0));
    base_req.stage3_gapopen = cli.has("-stage3_gapopen")
        ? static_cast<int16_t>(cli.get_int("-stage3_gapopen", 0))
        : INT16_MIN;
    base_req.stage3_gapext = cli.has("-stage3_gapext")
        ? static_cast<int16_t>(cli.get_int("-stage3_gapext", 0))
        : INT16_MIN;
    {
        double min_pident = cli.get_double("-stage3_min_pident", 0.0);
        base_req.stage3_min_pident_x100 = static_cast<uint16_t>(min_pident * 100.0);
    }
    base_req.stage3_min_nident = static_cast<uint32_t>(cli.get_int("-stage3_min_nident", 0));

    // Context
    {
        std::string context_str = cli.get_string("-context", "0");
        if (context_str.find('.') != std::string::npos) {
            double ratio = std::stod(context_str);
            base_req.context_frac_x10000 = static_cast<uint16_t>(ratio * 10000.0);
        } else {
            base_req.context_abs = static_cast<uint32_t>(std::stoi(context_str));
        }
    }

    // Validate output format compatibility
    {
        std::string err;
        if (!validate_output_format(outfmt, base_req.mode,
                                    base_req.stage3_traceback != 0,
                                    output_path, err)) {
            std::fprintf(stderr, "%s\n", err.c_str());
            return 1;
        }
    }

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

    // Pre-flight validation
    InfoResponse server_info;
    if (!execute_info(cli, has_http, server_info, logger
#ifdef IKAFSSN_ENABLE_HTTP
                      , auth
#endif
                      )) {
        return 1;
    }
    {
        std::string err = validate_info(server_info, base_req.db,
                                        base_req.k, base_req.mode, true);
        if (!err.empty()) {
            std::fprintf(stderr, "%s\n", err.c_str());
            return 1;
        }
    }
    logger.debug("Pre-flight validation passed");

    // Build query_id -> sequence lookup for retry
    std::unordered_map<std::string, std::string> query_map;
    for (const auto& q : queries) {
        query_map[q.id] = q.sequence;
    }

    // Determine batch size from server info
    int batch_size = static_cast<int>(queries.size());
    if (server_info.max_seqs_per_req > 0)
        batch_size = std::min(batch_size, static_cast<int>(server_info.max_seqs_per_req));
    if (server_info.max_queue_size > 0) {
        int available = server_info.max_queue_size - server_info.queue_depth;
        if (available > 0)
            batch_size = std::min(batch_size, available);
    }
    if (batch_size <= 0) batch_size = 1;
    logger.debug("Batch size: %d (queries=%zu, max_seqs_per_req=%d, available=%d/%d)",
                 batch_size, queries.size(), server_info.max_seqs_per_req,
                 server_info.max_queue_size - server_info.queue_depth,
                 server_info.max_queue_size);

    // Accumulate results across batches and retry attempts
    std::vector<OutputHit> all_hits;
    bool has_skipped = false;
    uint8_t resp_mode = 0;
    uint8_t resp_stage1_score = 0;
    bool resp_stage3_traceback = false;
    bool first_response = true;

    // Retry schedule: 30s, 60s, 120s, 120s, ...
    static constexpr int retry_delays[] = {30, 60, 120};
    static constexpr int num_retry_delays =
        static_cast<int>(sizeof(retry_delays) / sizeof(retry_delays[0]));

    // Helper lambda to collect results from a response
    auto collect_results = [&](const SearchResponse& resp) {
        for (const auto& qr : resp.results) {
            if (qr.skipped != 0) {
                has_skipped = true;
                std::fprintf(stderr, "Warning: query '%s' was skipped (degenerate bases)\n",
                             qr.qseqid.c_str());
                continue;
            }
            if (qr.warnings & kWarnMultiDegen) {
                std::fprintf(stderr,
                    "Warning: query '%s' contains k-mers with 2 or more degenerate bases; "
                    "those k-mers are ignored and not used in the search\n",
                    qr.qseqid.c_str());
            }
            for (const auto& hit : qr.hits) {
                OutputHit oh;
                oh.qseqid = qr.qseqid;
                oh.sseqid = hit.sseqid;
                oh.sstrand = (hit.sstrand == 0) ? '+' : '-';
                oh.qstart = hit.qstart;
                oh.qend = hit.qend;
                oh.sstart = hit.sstart;
                oh.send = hit.send;
                oh.chainscore = hit.chainscore;
                oh.coverscore = hit.coverscore;
                oh.matchscore = hit.matchscore;
                oh.volume = hit.volume;
                oh.qlen = hit.qlen;
                oh.slen = hit.slen;
                oh.alnscore = hit.alnscore;
                oh.nident = hit.nident;
                oh.mismatch = hit.mismatch;
                oh.pident = static_cast<double>(hit.pident_x100) / 100.0;
                oh.cigar = hit.cigar;
                oh.qseq = hit.qseq;
                oh.sseq = hit.sseq;
                all_hits.push_back(std::move(oh));
            }
        }
    };

    // Process queries in batches
    size_t sent = 0;
    while (sent < queries.size()) {
        size_t batch_end = std::min(sent + static_cast<size_t>(batch_size), queries.size());

        // Build request for this batch
        SearchRequest req = base_req;
        for (size_t i = sent; i < batch_end; i++) {
            req.queries.push_back({queries[i].id, queries[i].sequence});
        }
        sent = batch_end;

        logger.info("Sending batch of %zu queries (%zu/%zu)",
                     req.queries.size(), sent, queries.size());

        // Retry loop for this batch (rejected = server busy)
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
                         resp.k, resp.results.size(), resp.rejected_qseqids.size());

            // Save mode/score_type from first successful response
            if (first_response) {
                resp_mode = resp.mode;
                resp_stage1_score = resp.stage1_score;
                resp_stage3_traceback = (resp.stage3_traceback != 0);
                first_response = false;
            }

            collect_results(resp);

            // If no rejected queries, this batch is done
            if (resp.rejected_qseqids.empty()) break;

            // Sleep before retry
            int delay_idx = std::min(attempt, num_retry_delays - 1);
            int delay = retry_delays[delay_idx];
            logger.info("%zu queries rejected, retrying in %d seconds...",
                        resp.rejected_qseqids.size(), delay);
            std::this_thread::sleep_for(std::chrono::seconds(delay));

            // Rebuild request with only rejected query IDs
            req = base_req;
            for (const auto& qid : resp.rejected_qseqids) {
                auto it = query_map.find(qid);
                if (it != query_map.end()) {
                    req.queries.push_back({qid, it->second});
                }
            }

            if (req.queries.empty()) break; // all resolved
        }
    }

    // Write output using mode/stage1_score from response
    if (!write_all_results(output_path, all_hits, outfmt,
                           resp_mode, resp_stage1_score,
                           resp_stage3_traceback)) {
        return 1;
    }

    logger.info("Done. %zu hit(s) reported.", all_hits.size());
    return has_skipped ? 2 : 0;
}
