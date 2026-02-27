#include "ikafssnclient/socket_client.hpp"
#ifdef IKAFSSN_ENABLE_HTTP
#include "ikafssnclient/http_client.hpp"
#endif
#include "ikafssnclient/checkpoint.hpp"
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
#include <iostream>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>

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
        "  -ix <name>               Target database name on server\n"
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

// Convert SearchResponse results into OutputHit vector
static void collect_results(const SearchResponse& resp,
                             std::vector<OutputHit>& hits,
                             bool& has_skipped) {
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
            hits.push_back(std::move(oh));
        }
    }
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

    if (!cli.has("-ix")) {
        std::fprintf(stderr, "Error: -ix is required\n");
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
    std::string ix_name = cli.get_string("-ix");

    // Output format
    OutputFormat outfmt;
    {
        std::string err;
        if (!parse_output_format(cli.get_string("-outfmt", "tab"), outfmt, err)) {
            std::fprintf(stderr, "%s\n", err.c_str());
            return 1;
        }
    }

    // Read query FASTA (with stdin buffering for checkpoint)
    std::string stdin_content;
    std::vector<FastaRecord> queries;
    if (query_path == "-") {
        std::ostringstream oss;
        oss << std::cin.rdbuf();
        stdin_content = oss.str();
        std::istringstream iss(stdin_content);
        queries = read_fasta_stream(iss);
    } else {
        queries = read_fasta(query_path);
    }
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
    base_req.db = ix_name;

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
    std::string seqidlist_path;
    std::string neg_seqidlist_path;
    if (cli.has("-seqidlist")) {
        base_req.seqidlist_mode = SeqidlistMode::kInclude;
        seqidlist_path = cli.get_string("-seqidlist");
        base_req.seqids = read_seqidlist(seqidlist_path);
        logger.info("Loaded %zu accessions (include mode)", base_req.seqids.size());
    } else if (cli.has("-negative_seqidlist")) {
        base_req.seqidlist_mode = SeqidlistMode::kExclude;
        neg_seqidlist_path = cli.get_string("-negative_seqidlist");
        base_req.seqids = read_seqidlist(neg_seqidlist_path);
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

    // Resolve k value from server info
    uint8_t resolved_k = base_req.k;
    if (resolved_k == 0) {
        for (const auto& db : server_info.databases) {
            if (db.name == base_req.db) {
                resolved_k = db.default_k;
                break;
            }
        }
    }

    // Compute SHA256s for checkpoint validation
    std::string input_sha256;
    if (query_path == "-") {
        input_sha256 = sha256_string(stdin_content);
    } else {
        input_sha256 = sha256_file(query_path);
    }
    std::string seqidlist_sha256;
    std::string neg_seqidlist_sha256;
    if (!seqidlist_path.empty())
        seqidlist_sha256 = sha256_file(seqidlist_path);
    if (!neg_seqidlist_path.empty())
        neg_seqidlist_sha256 = sha256_file(neg_seqidlist_path);

    // Build options text and checkpoint
    DbStats db_stats = resolve_db_stats(server_info, base_req.db, resolved_k);
    std::string options_text = build_options_text(
        base_req, db_stats, resolved_k, outfmt,
        seqidlist_sha256, neg_seqidlist_sha256);

    Checkpoint::Config ckpt_cfg;
    ckpt_cfg.output_path = output_path;
    ckpt_cfg.input_path = query_path;
    ckpt_cfg.ix_name = ix_name;
    ckpt_cfg.resolved_k = resolved_k;
    ckpt_cfg.outfmt = outfmt;

    Checkpoint ckpt(ckpt_cfg, logger);

    // Acquire lock
    LockGuard lock;
    if (!ckpt.acquire_lock(lock)) {
        // Lock dir might not exist yet if temp_dir doesn't exist
        // Try to create temp dir first, then lock
    }

    // Resume or initialize checkpoint
    std::unordered_set<std::string> completed_seqids;
    int next_batch_num = 0;

    if (ckpt.exists()) {
        // Acquire lock first if we haven't already
        if (!lock.locked()) {
            if (!ckpt.acquire_lock(lock)) {
                return 1;
            }
        }
        if (!ckpt.resume(options_text, input_sha256,
                          completed_seqids, next_batch_num)) {
            logger.info("Checkpoint validation failed, starting fresh");
            lock.release();
            ckpt.cleanup();
            if (!ckpt.initialize(options_text, input_sha256, stdin_content)) {
                return 1;
            }
            if (!ckpt.acquire_lock(lock)) {
                return 1;
            }
        } else {
            logger.info("Resumed from checkpoint: %zu queries already completed",
                        completed_seqids.size());
        }
    } else {
        if (!ckpt.initialize(options_text, input_sha256, stdin_content)) {
            return 1;
        }
        if (!lock.locked()) {
            if (!ckpt.acquire_lock(lock)) {
                return 1;
            }
        }
    }

    // Build remaining queries (filter out completed)
    std::vector<FastaRecord> remaining;
    for (const auto& q : queries) {
        if (completed_seqids.find(q.id) == completed_seqids.end()) {
            remaining.push_back(q);
        }
    }
    logger.info("%zu remaining queries to process", remaining.size());

    // Build query_id -> sequence lookup for retry
    std::unordered_map<std::string, std::string> query_map;
    for (const auto& q : queries) {
        query_map[q.id] = q.sequence;
    }

    // Determine batch size from server info
    int batch_size = static_cast<int>(remaining.size());
    if (batch_size == 0) batch_size = 1;  // avoid div-by-zero
    if (server_info.max_seqs_per_req > 0)
        batch_size = std::min(batch_size, static_cast<int>(server_info.max_seqs_per_req));
    if (server_info.max_queue_size > 0) {
        int available = server_info.max_queue_size - server_info.queue_depth;
        if (available > 0)
            batch_size = std::min(batch_size, available);
    }
    if (batch_size <= 0) batch_size = 1;
    logger.debug("Batch size: %d (remaining=%zu, max_seqs_per_req=%d)",
                 batch_size, remaining.size(), server_info.max_seqs_per_req);

    bool has_skipped = false;
    uint8_t resp_mode = 0;
    uint8_t resp_stage1_score = 0;
    bool resp_stage3_traceback = false;
    bool first_response = true;

    // Try to load response meta from checkpoint (for resume with all queries done)
    if (!remaining.empty() || !ckpt.read_response_meta(resp_mode, resp_stage1_score,
                                                        resp_stage3_traceback)) {
        // Will be populated from first response
    } else {
        first_response = false;
    }

    // Retry schedule: 30s, 60s, 120s, 120s, ...
    static constexpr int retry_delays[] = {30, 60, 120};
    static constexpr int num_retry_delays =
        static_cast<int>(sizeof(retry_delays) / sizeof(retry_delays[0]));

    int batch_num = next_batch_num;

    // Process remaining queries in batches
    size_t sent = 0;
    while (sent < remaining.size()) {
        size_t batch_end = std::min(sent + static_cast<size_t>(batch_size),
                                     remaining.size());

        // Build batch seqid list
        std::vector<std::string> batch_seqids;
        for (size_t i = sent; i < batch_end; i++) {
            batch_seqids.push_back(remaining[i].id);
        }

        // Write batch seqids before sending
        ckpt.write_batch_seqids(batch_num, batch_seqids);

        // Build request for this batch
        SearchRequest req = base_req;
        for (size_t i = sent; i < batch_end; i++) {
            req.queries.push_back({remaining[i].id, remaining[i].sequence});
        }
        sent = batch_end;

        logger.info("Sending batch %d: %zu queries (%zu/%zu)",
                     batch_num, req.queries.size(), sent, remaining.size());

        // Retry loop for this batch (rejected = server busy)
        for (int attempt = 0; ; attempt++) {
            SearchResponse resp;
            if (!execute_search(cli, has_http, req, resp, logger
#ifdef IKAFSSN_ENABLE_HTTP
                    , auth
#endif
                    )) {
                lock.release();
                return 1;
            }

            if (resp.status != 0) {
                std::fprintf(stderr, "Error: server returned status %d\n", resp.status);
                lock.release();
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
                ckpt.write_response_meta(resp_mode, resp_stage1_score,
                                          resp_stage3_traceback);
            }

            // Collect accepted results
            std::vector<OutputHit> batch_hits;
            collect_results(resp, batch_hits, has_skipped);

            if (resp.rejected_qseqids.empty()) {
                // All accepted - save batch results
                ckpt.write_batch_results(batch_num, batch_hits,
                                          resp_mode, resp_stage1_score,
                                          resp_stage3_traceback);
                batch_num++;
                break;
            }

            // Partial reject: save accepted results with updated seqid list
            std::unordered_set<std::string> rejected_set(
                resp.rejected_qseqids.begin(), resp.rejected_qseqids.end());

            // Rewrite batch seqids to only include accepted
            std::vector<std::string> accepted_seqids;
            for (const auto& id : batch_seqids) {
                if (rejected_set.find(id) == rejected_set.end()) {
                    accepted_seqids.push_back(id);
                }
            }
            ckpt.write_batch_seqids(batch_num, accepted_seqids);
            ckpt.write_batch_results(batch_num, batch_hits,
                                      resp_mode, resp_stage1_score,
                                      resp_stage3_traceback);
            batch_num++;

            // Sleep before retry
            int delay_idx = std::min(attempt, num_retry_delays - 1);
            int delay = retry_delays[delay_idx];
            logger.info("%zu queries rejected, retrying in %d seconds...",
                        resp.rejected_qseqids.size(), delay);
            std::this_thread::sleep_for(std::chrono::seconds(delay));

            // Rebuild request with only rejected query IDs as new batch
            batch_seqids.clear();
            req = base_req;
            for (const auto& qid : resp.rejected_qseqids) {
                auto it = query_map.find(qid);
                if (it != query_map.end()) {
                    req.queries.push_back({qid, it->second});
                    batch_seqids.push_back(qid);
                }
            }

            if (req.queries.empty()) break; // all resolved

            // Write seqids for the new retry batch
            ckpt.write_batch_seqids(batch_num, batch_seqids);
        }
    }

    // If all queries were already completed (resume case), load meta
    if (remaining.empty() && first_response) {
        if (!ckpt.read_response_meta(resp_mode, resp_stage1_score,
                                      resp_stage3_traceback)) {
            logger.error("No response metadata found in checkpoint");
            lock.release();
            return 1;
        }
    }

    // Merge all batch results
    if (!ckpt.merge_results(output_path, resp_mode, resp_stage1_score,
                             resp_stage3_traceback)) {
        logger.error("Failed to merge results");
        lock.release();
        return 1;
    }

    // Cleanup
    lock.release();
    ckpt.cleanup();

    logger.info("Done.");
    return has_skipped ? 2 : 0;
}
