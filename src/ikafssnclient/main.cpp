#include "ikafssnclient/socket_client.hpp"
#ifdef IKAFSSN_ENABLE_HTTP
#include "ikafssnclient/http_client.hpp"
#include "ikafssnclient/async_http_client.hpp"
#include "ikafssnclient/job_dir.hpp"
#include "ikafssnclient/poll_loop.hpp"
#include "ikafssnclient/failed_writer.hpp"
#endif
#include "ikafssnclient/checkpoint.hpp"
#include "core/spaced_seed.hpp"
#include "core/version.hpp"
#include "protocol/info_format.hpp"
#include "protocol/messages.hpp"
#include "protocol/serializer.hpp"
#include "util/cli_validators.hpp"
#include "util/common_init.hpp"
#include "util/simd_dispatch.hpp"
#include "io/compressed_stream.hpp"
#include "io/fasta_reader.hpp"
#include "io/primer_query.hpp"
#include "io/seqidlist_reader.hpp"
#include "io/result_writer.hpp"
#include "io/sam_writer.hpp"
#include "util/cli_parser.hpp"
#include "util/socket_utils.hpp"
#include "util/logger.hpp"

#include <chrono>
#include <climits>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>

using namespace ikafssn;

#ifdef IKAFSSN_ENABLE_HTTP
static void install_sigint_handler() {
    static bool installed = false;
    if (installed) return;
    installed = true;
    std::signal(SIGINT, [](int) {
        PollLoop::sigint_flag().store(true);
    });
}
#endif

static void print_usage(const char* prog) {
    print_version_header("ikafssnclient");
    std::fprintf(stderr,
        "Usage: %s [options]\n"
        "\n"
        "Connection (one required for search):\n"
        "  -socket <path>           UNIX domain socket path\n"
        "  -tcp <host>:<port>       TCP server address\n"
#ifdef IKAFSSN_ENABLE_HTTP
        "  -http <url>              ikafssnhttpd URL (e.g., http://example.com:8080)\n"
        "                           HTTP mode uses async REST polling.\n"
#endif
        "\n"
        "Required:\n"
        "  -query <path>            Query FASTA file (- for stdin)\n"
        "  -ix <name>               Target database name on server\n"
        "\n"
        "Filtering:\n"
        "  -min_query_length <int>  Minimum query length; shorter queries are skipped\n"
        "                           (default: 64; must be >= server's min_seq_length)\n"
        "\n"
#ifdef IKAFSSN_ENABLE_HTTP
        "Async REST job management (requires -http):\n"
        "  -submit_only             Submit, print group_id, and exit\n"
        "  -jobs                    List all locally-tracked job groups\n"
        "  -job_detail <id>         Show jobs in a group, or detail of a single job\n"
        "  -resume <id>             Resume polling for an existing group or job\n"
#endif
        "(See `ikafssnclient -h` for the full list of search options.)\n",
        prog);
}

#ifdef IKAFSSN_ENABLE_HTTP
// Inline subcommand handlers for the async REST stack.

static int cmd_jobs() {
    auto root = default_jobs_root();
    auto groups = list_groups(root);
    if (groups.empty()) {
        std::printf("No locally-tracked job groups under %s\n", root.c_str());
        return 0;
    }
    std::printf("%-40s  %-20s  %s\n", "GROUP_ID", "SUBMITTED_AT", "JOBS (q/r/d/f)");
    for (const auto& g : groups) {
        GroupMeta gm;
        std::string err;
        if (!read_group_meta(root, g, gm, err)) continue;
        char tbuf[32];
        time_t t = static_cast<time_t>(gm.submitted_at);
        std::strftime(tbuf, sizeof(tbuf), "%Y-%m-%dT%H:%M:%S",
                      std::localtime(&t));
        std::printf("%-40s  %-20s  %d/%d/%d/%d\n",
                    gm.group_id.c_str(), tbuf,
                    gm.cnt_queued, gm.cnt_running, gm.cnt_done, gm.cnt_failed);
    }
    return 0;
}

static int cmd_jobdetail(const std::string& id) {
    auto root = default_jobs_root();
    bool is_group = false;
    std::string gid, jid;
    if (!resolve_id(root, id, is_group, gid, jid)) {
        std::fprintf(stderr, "Error: id not found: %s\n", id.c_str());
        return 1;
    }
    if (is_group) {
        GroupMeta gm;
        std::string err;
        if (!read_group_meta(root, gid, gm, err)) {
            std::fprintf(stderr, "Error: %s\n", err.c_str());
            return 1;
        }
        std::printf("group_id=%s db=%s url=%s\n",
                    gm.group_id.c_str(), gm.db.c_str(), gm.httpd_url.c_str());
        std::printf("query=%s output=%s output_format=%s\n",
                    gm.query_file_path_abs.c_str(),
                    gm.output_path.c_str(), gm.output_format.c_str());
        std::printf("%-40s  %-10s  %-8s\n", "JOB_ID", "STATUS", "ATTEMPTS");
        for (const auto& job_id : gm.job_ids) {
            JobMeta jm;
            std::string e;
            if (!read_job_meta(root, gid, job_id, jm, e)) continue;
            std::printf("%-40s  %-10s  %-8d\n",
                        jm.job_id.c_str(), jm.status.c_str(), jm.attempts);
        }
    } else {
        JobMeta jm;
        std::string err;
        if (!read_job_meta(root, gid, jid, jm, err)) {
            std::fprintf(stderr, "Error: %s\n", err.c_str());
            return 1;
        }
        std::printf("job_id=%s\ngroup_id=%s\nstatus=%s\nattempts=%d\n",
                    jm.job_id.c_str(), jm.group_id.c_str(),
                    jm.status.c_str(), jm.attempts);
        std::printf("submitted_at=%lld completed_at=%lld\n",
                    (long long)jm.submitted_at, (long long)jm.completed_at);
        if (!jm.error_message.empty())
            std::printf("error_message=%s\n", jm.error_message.c_str());
        if (!jm.fail_reason.empty())
            std::printf("fail_reason=%s\n", jm.fail_reason.c_str());
    }
    return 0;
}

// After polling completes, walk every job in the group, deserialise its
// cached result blob, and merge into the user's output file.  Failed
// jobs are turned into kFailHttpJob OutputHit rows via failed_writer.
static int finalize_group(const std::string& root,
                          const GroupMeta& gm,
                          Logger& logger) {
    std::vector<OutputHit> all_hits;
    uint8_t resp_mode = 0;
    uint8_t resp_stage1 = 0;
    bool resp_traceback = false;
    bool first_resp = true;

    for (const auto& job_id : gm.job_ids) {
        JobMeta jm;
        std::string err;
        if (!read_job_meta(root, gm.group_id, job_id, jm, err)) {
            logger.error("finalize: missing job meta %s: %s",
                         job_id.c_str(), err.c_str());
            continue;
        }

        if (jm.status == "done") {
            std::vector<uint8_t> blob;
            if (!read_job_result(root, gm.group_id, job_id, blob, err)) {
                logger.error("finalize: read result %s: %s",
                             job_id.c_str(), err.c_str());
                continue;
            }
            SearchResponse resp;
            if (!deserialize(blob, resp)) {
                logger.error("finalize: deserialize %s failed",
                             job_id.c_str());
                continue;
            }
            if (first_resp) {
                resp_mode = resp.mode;
                resp_stage1 = resp.stage1_score;
                resp_traceback = (resp.stage3_traceback != 0);
                first_resp = false;
            }
            for (const auto& qr : resp.results) {
                if (qr.skip_reason != 0) {
                    OutputHit oh;
                    oh.qseqid = qr.qseqid;
                    oh.qlen = qr.qlen;
                    oh.skip_reason = qr.skip_reason;
                    oh.skip_detail = qr.skip_detail;
                    oh.sstrand = '*';
                    all_hits.push_back(std::move(oh));
                    continue;
                }
                for (const auto& hit : qr.hits) {
                    OutputHit oh;
                    oh.qseqid = qr.qseqid;
                    oh.sseqid = hit.sseqid;
                    oh.sstrand = (hit.sstrand == 0) ? '+' : '-';
                    oh.qstart = hit.qstart; oh.qend = hit.qend;
                    oh.sstart = hit.sstart; oh.send = hit.send;
                    oh.chainscore = hit.chainscore;
                    oh.coverscore = hit.coverscore;
                    oh.matchscore = hit.matchscore;
                    oh.volume = hit.volume;
                    oh.qlen = hit.qlen; oh.slen = hit.slen;
                    oh.alnscore = hit.alnscore;
                    oh.npositive = hit.npositive;
                    oh.nnegative = hit.nnegative;
                    oh.ppositive =
                        static_cast<double>(hit.ppositive_x100) / 100.0;
                    oh.cigar = hit.cigar;
                    oh.qseq = hit.qseq; oh.sseq = hit.sseq;
                    all_hits.push_back(std::move(oh));
                }
            }
        } else if (jm.status == "failed" || jm.status == "timeout") {
            std::vector<std::string> deflines;
            std::string derr;
            if (!read_job_deflines(root, gm.group_id, job_id, deflines, derr)) {
                logger.error("finalize: deflines %s: %s",
                             job_id.c_str(), derr.c_str());
                continue;
            }
            std::string reason = !jm.fail_reason.empty()
                ? jm.fail_reason
                : (jm.status == "timeout"
                    ? "timeout: exceeded retention_time"
                    : "backend_error: unknown");
            synth_failed_hits(deflines, reason, all_hits);
        }
    }

    OutputFormat fmt = OutputFormat::kTsv;
    if (gm.output_format == "json") fmt = OutputFormat::kJson;
    else if (gm.output_format == "sam") fmt = OutputFormat::kSam;
    else if (gm.output_format == "bam") fmt = OutputFormat::kBam;

    if (!write_all_results(gm.output_path, all_hits, fmt,
                            resp_mode, resp_stage1, resp_traceback,
                            gm.compression_level)) {
        logger.error("finalize: write_all_results failed");
        return 1;
    }
    // Drop result.bin.zst caches; metadata stays for forensic use.
    for (const auto& job_id : gm.job_ids) {
        std::string e;
        delete_job_result(root, gm.group_id, job_id, e);
    }
    return 0;
}

static int cmd_resume(const std::string& id, Logger& logger,
                      const HttpAuthConfig& auth) {
    install_sigint_handler();
    auto root = default_jobs_root();
    bool is_group = false;
    std::string gid, jid;
    if (!resolve_id(root, id, is_group, gid, jid)) {
        std::fprintf(stderr, "Error: id not found: %s\n", id.c_str());
        return 1;
    }
    GroupMeta gm;
    std::string err;
    if (!read_group_meta(root, gid, gm, err)) {
        std::fprintf(stderr, "Error: %s\n", err.c_str());
        return 1;
    }
    PollLoop loop(root, gid, auth, logger);
    if (!loop.run()) {
        if (PollLoop::sigint_flag().load()) {
            return 130;
        }
        // Other failure: fall through to finalize what we have.
    }
    return finalize_group(root, gm, logger);
}
#endif // IKAFSSN_ENABLE_HTTP

// Execute an info request via socket / TCP / HTTP. Returns true on success.
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
    } else {
        std::string tcp_addr = cli.get_string("-tcp");
        fd = tcp_connect(tcp_addr);
    }
    if (fd < 0) {
        std::fprintf(stderr, "Error: cannot connect for info\n");
        return false;
    }
    if (!socket_info(fd, resp)) {
        std::fprintf(stderr, "Error: info request failed\n");
        close_fd(fd);
        return false;
    }
    close_fd(fd);
    return true;
}

static bool execute_search(
    const CliParser& cli,
    const SearchRequest& req,
    SearchResponse& resp,
    const Logger& logger) {
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
        if (qr.skip_reason != 0) {
            has_skipped = true;
            std::fprintf(stderr, "Warning: query '%s' skipped: %s (%s)\n",
                         qr.qseqid.c_str(),
                         skip_reason_str(qr.skip_reason),
                         qr.skip_detail.c_str());
            OutputHit oh;
            oh.qseqid = qr.qseqid;
            oh.qlen = qr.qlen;
            oh.skip_reason = qr.skip_reason;
            oh.skip_detail = qr.skip_detail;
            oh.sstrand = '*';
            hits.push_back(std::move(oh));
            continue;
        }
        if (qr.warnings & kWarnMultiDegen) {
            std::fprintf(stderr,
                "Warning: query '%s' contains k-mers with 2 or more degenerate bases\n",
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
            oh.npositive = hit.npositive;
            oh.nnegative = hit.nnegative;
            oh.ppositive = static_cast<double>(hit.ppositive_x100) / 100.0;
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

    if (cli.has("-h") || cli.has("-help")) {
        print_usage(argv[0]);
        return 0;
    }

    bool has_http = false;
#ifdef IKAFSSN_ENABLE_HTTP
    has_http = cli.has("-http");
#endif

#ifdef IKAFSSN_ENABLE_HTTP
    Logger early_logger = make_logger(cli);
    HttpAuthConfig auth;
    if (cli.has("-user") && cli.has("-http_user")) {
        std::fprintf(stderr, "Error: -user and -http_user are mutually exclusive\n");
        return 1;
    }
    if (cli.has("-user")) auth.userpwd = cli.get_string("-user");
    else if (cli.has("-http_user")) {
        std::string user = cli.get_string("-http_user");
        std::string pass = cli.get_string("-http_password", "");
        auth.userpwd = user + ":" + pass;
    }
    if (cli.has("-netrc_file")) auth.netrc_file = cli.get_string("-netrc_file");

    if (cli.has("-jobs")) {
        return cmd_jobs();
    }
    if (cli.has("-job_detail")) {
        return cmd_jobdetail(cli.get_string("-job_detail"));
    }
    if (cli.has("-resume")) {
        return cmd_resume(cli.get_string("-resume"), early_logger, auth);
    }
    if (cli.has("-submit_only") && !has_http) {
        std::fprintf(stderr, "Error: -submit_only requires -http\n");
        return 1;
    }
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

    bool has_query = cli.has("-query");
    bool has_primer = cli.has("-primer");
    if (!has_query && !has_primer) {
        std::fprintf(stderr, "Error: -query or -primer is required\n");
        return 1;
    }
    if (has_query && has_primer) {
        std::fprintf(stderr, "Error: -query and -primer are mutually exclusive\n");
        return 1;
    }
    {
        std::string err;
        if (!validate_primer_mode_options(cli, err)) {
            std::fprintf(stderr, "%s\n", err.c_str());
            return 1;
        }
    }
    if (!cli.has("-ix")) {
        std::fprintf(stderr, "Error: -ix is required\n");
        return 1;
    }
    if (cli.has("-seqidlist") && cli.has("-negative_seqidlist")) {
        std::fprintf(stderr, "Error: -seqidlist and -negative_seqidlist are mutually exclusive\n");
        return 1;
    }

    Logger logger = make_logger(cli);
    init_simd_dispatch(&logger);

    std::string query_path = has_query ? cli.get_string("-query") : "";
    std::string primer_path = has_primer ? cli.get_string("-primer") : "";
    std::string input_path = has_primer ? primer_path : query_path;
    std::string output_path = cli.get_string("-o");
    std::string ix_name = cli.get_string("-ix");

    uint32_t insert_length = 0;
    double stage1_primer_score = 0.5;
    int stage2_primer_score_add = 1;
    if (has_primer) {
        insert_length = static_cast<uint32_t>(cli.get_int("-insert_length", 0));
        stage1_primer_score = cli.get_double("-stage1_primer_score", 0.5);
        stage2_primer_score_add = cli.get_int("-stage2_primer_score_add", 1);
        if (stage1_primer_score > 1.0 && stage1_primer_score < 2.0) {
            std::fprintf(stderr, "Error: -stage1_primer_score must be 0<v<=1 or v>=2\n");
            return 1;
        }
    }

    OutputFormat outfmt;
    {
        std::string err;
        if (!parse_output_format(cli.get_string("-output_format", "tsv"), outfmt, err)) {
            std::fprintf(stderr, "%s\n", err.c_str());
            return 1;
        }
    }

    std::string stdin_content;
    std::vector<FastaRecord> queries;
    if (input_path == "-") {
        std::ostringstream oss;
        oss << std::cin.rdbuf();
        stdin_content = oss.str();
        std::istringstream iss(stdin_content);
        queries = read_fasta_stream(iss);
    } else {
        queries = read_fasta(input_path);
    }
    if (queries.empty()) {
        std::fprintf(stderr, "Error: no %s sequences found\n",
                     has_primer ? "primer" : "query");
        return 1;
    }

    // -min_query_length default 64.  The integrity check against the
    // server's min_seq_length runs after execute_info().
    uint32_t min_query_length = 64;
    {
        int v = cli.get_int("-min_query_length", 64);
        if (v <= 0) {
            std::fprintf(stderr,
                "Error: -min_query_length must be a positive integer (got %d)\n", v);
            return 1;
        }
        min_query_length = static_cast<uint32_t>(v);
    }

    SearchRequest base_req;
    base_req.k = static_cast<uint8_t>(cli.get_int("-k", 0));
    if (cli.has("-stage2_min_score")) {
        base_req.stage2_min_score = static_cast<uint16_t>(cli.get_int("-stage2_min_score", 0));
        base_req.has_stage2_min_score = 1;
    }
    base_req.stage2_max_gap = static_cast<uint16_t>(cli.get_int("-stage2_max_gap", 0));
    base_req.stage2_max_lookback = static_cast<uint16_t>(cli.get_int("-stage2_max_lookback", 0));
    base_req.stage2_max_nhit_per_subject = static_cast<uint16_t>(cli.get_int("-stage2_max_nhit_per_subject", 0));
    base_req.stage2_min_nhit_diag = static_cast<uint8_t>(cli.get_int("-stage2_min_nhit_diag", 0));
    base_req.stage1_topn = static_cast<uint16_t>(cli.get_int("-stage1_topn", 0));
    {
        double min_s1 = cli.get_double("-stage1_min_score", 0.0);
        if (min_s1 > 0 && min_s1 < 1.0) {
            base_req.stage1_min_score_frac_x10000 = static_cast<uint16_t>(min_s1 * 10000.0);
        } else {
            base_req.stage1_min_score = static_cast<uint16_t>(min_s1);
        }
    }
    base_req.nresult = static_cast<uint16_t>(cli.get_int("-nresult", 0));
    base_req.mode = static_cast<uint8_t>(cli.get_int("-mode", 0));
    base_req.stage1_score = static_cast<uint8_t>(cli.get_int("-stage1_score", 0));
    base_req.accept_qdegen = static_cast<uint8_t>(cli.get_int("-accept_qdegen", 1));
    base_req.strand = static_cast<int8_t>(cli.get_int("-strand", 0));
    base_req.db = ix_name;
    {
        std::string err;
        if (!parse_max_degen_expand(cli, 0, base_req.max_degen_expand, err)) {
            std::fprintf(stderr, "%s\n", err.c_str());
            return 1;
        }
    }
    {
        std::string err;
        if (!parse_spaced_seed_t(cli, base_req.t, err)) {
            std::fprintf(stderr, "%s\n", err.c_str());
            return 1;
        }
    }
    if (cli.has("-template_type")) {
        TemplateType tt;
        std::string err;
        if (!parse_template_type_cli(cli, TemplateType::kBoth, false, tt, err)) {
            std::fprintf(stderr, "%s\n", err.c_str());
            return 1;
        }
        base_req.template_type = static_cast<uint8_t>(tt);
    }
    base_req.stage3_traceback = static_cast<uint8_t>(cli.get_int("-stage3_traceback", 0));
    base_req.stage3_gapopen = cli.has("-stage3_gapopen")
        ? static_cast<int16_t>(cli.get_int("-stage3_gapopen", 0)) : INT16_MIN;
    base_req.stage3_gapext = cli.has("-stage3_gapext")
        ? static_cast<int16_t>(cli.get_int("-stage3_gapext", 0)) : INT16_MIN;
    {
        double min_ppositive = cli.get_double("-stage3_min_ppositive", 0.0);
        base_req.stage3_min_ppositive_x100 = static_cast<uint16_t>(min_ppositive * 100.0);
    }
    base_req.stage3_min_npositive = static_cast<uint32_t>(cli.get_int("-stage3_min_npositive", 0));
    {
        std::string sm;
        std::string err;
        if (!parse_score_matrix(cli, sm, err)) {
            std::fprintf(stderr, "%s\n", err.c_str());
            return 1;
        }
        base_req.score_matrix = score_matrix_code(sm);
    }
    {
        std::string context_str = cli.get_string("-context_extend", "2.0");
        if (context_str.find('.') != std::string::npos) {
            double ratio = std::stod(context_str);
            base_req.context_frac_x10000 = static_cast<uint16_t>(ratio * 10000.0);
        } else {
            base_req.context_abs = static_cast<uint32_t>(std::stoi(context_str));
        }
    }
    {
        std::string err;
        if (!validate_output_format(outfmt, base_req.mode,
                                    base_req.stage3_traceback != 0,
                                    output_path, err)) {
            std::fprintf(stderr, "%s\n", err.c_str());
            return 1;
        }
    }

    int compression_level = cli.has("-compression_level")
        ? cli.get_int("-compression_level", kCompressionLevelDefault)
        : kCompressionLevelDefault;
    {
        std::string err;
        if (!validate_compression_level(
                detect_format_from_extension(output_path),
                compression_level, err)) {
            std::fprintf(stderr, "%s\n", err.c_str());
            return 1;
        }
    }

    std::string seqidlist_path, neg_seqidlist_path;
    if (cli.has("-seqidlist")) {
        base_req.seqidlist_mode = SeqidlistMode::kInclude;
        seqidlist_path = cli.get_string("-seqidlist");
        base_req.seqids = read_seqidlist(seqidlist_path);
    } else if (cli.has("-negative_seqidlist")) {
        base_req.seqidlist_mode = SeqidlistMode::kExclude;
        neg_seqidlist_path = cli.get_string("-negative_seqidlist");
        base_req.seqids = read_seqidlist(neg_seqidlist_path);
    }

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
                                        base_req.k, base_req.mode, true,
                                        base_req.t, base_req.template_type);
        if (!err.empty()) {
            std::fprintf(stderr, "%s\n", err.c_str());
            return 1;
        }
    }

    // Pre-filter queries by -min_query_length and the server-reported
    // overlap_length.  The server enforces the same checks but rejecting
    // here saves a round trip and surfaces a clear local error.
    {
        uint32_t srv_overlap_length = 0;
        for (const auto& db : server_info.databases) {
            if (db.name == base_req.db) {
                srv_overlap_length = db.overlap_length;
                break;
            }
        }

        std::vector<FastaRecord> kept;
        kept.reserve(queries.size());
        uint32_t skipped_short = 0;
        uint32_t skipped_long  = 0;
        for (auto& q : queries) {
            if (q.sequence.size() < min_query_length) {
                std::fprintf(stderr,
                    "Warning: query '%s' (length=%zu) is shorter than -min_query_length=%u, "
                    "skipped\n",
                    q.id.c_str(), q.sequence.size(),
                    static_cast<unsigned>(min_query_length));
                ++skipped_short;
                continue;
            }
            if (srv_overlap_length > 0 &&
                q.sequence.size() > srv_overlap_length) {
                std::fprintf(stderr,
                    "Warning: query '%s' (length=%zu) is longer than the index's "
                    "overlap_length=%u for database '%s', skipped\n",
                    q.id.c_str(), q.sequence.size(),
                    static_cast<unsigned>(srv_overlap_length),
                    base_req.db.c_str());
                ++skipped_long;
                continue;
            }
            kept.push_back(std::move(q));
        }
        if (skipped_short > 0) {
            logger.info("Skipped %u queries shorter than min_query_length=%u",
                        skipped_short, static_cast<unsigned>(min_query_length));
        }
        if (skipped_long > 0) {
            logger.info("Skipped %u queries longer than overlap_length=%u",
                        skipped_long, static_cast<unsigned>(srv_overlap_length));
        }
        queries = std::move(kept);
        if (queries.empty()) {
            std::fprintf(stderr,
                "Error: all input queries were rejected by -min_query_length=%u "
                "or the index's overlap_length=%u\n",
                static_cast<unsigned>(min_query_length),
                static_cast<unsigned>(srv_overlap_length));
            return 1;
        }
    }

    uint8_t resolved_k = base_req.k;
    if (resolved_k == 0) {
        for (const auto& db : server_info.databases) {
            if (db.name == base_req.db) {
                resolved_k = db.default_k;
                break;
            }
        }
    }

    if (has_primer) {
        std::vector<uint32_t> seed_masks;
        if (base_req.t > 0) {
            if (base_req.template_type == 3) {
                auto cod = get_seed_masks(resolved_k, base_req.t, TemplateType::kCoding);
                auto opt = get_seed_masks(resolved_k, base_req.t, TemplateType::kOptimal);
                seed_masks.insert(seed_masks.end(), cod.begin(), cod.end());
                seed_masks.insert(seed_masks.end(), opt.begin(), opt.end());
            } else {
                seed_masks = get_seed_masks(resolved_k, base_req.t,
                                 static_cast<TemplateType>(base_req.template_type));
            }
        }
        PrimerConfig pcfg;
        pcfg.insert_length = insert_length;
        pcfg.k = resolved_k;
        pcfg.t = base_req.t;
        pcfg.masks = base_req.t > 0 ? &seed_masks : nullptr;
        std::vector<PrimerPair> primer_pairs;
        std::string primer_err = parse_primer_pairs(queries, pcfg, primer_pairs);
        if (!primer_err.empty()) {
            std::fprintf(stderr, "%s\n", primer_err.c_str());
            return 1;
        }
        queries.clear();
        uint32_t min_stage2_score = UINT32_MAX;
        for (const auto& pp : primer_pairs) {
            FastaRecord qr;
            qr.id = pp.query_id;
            qr.sequence = pp.query_seq;
            queries.push_back(std::move(qr));
            uint32_t total_pos = pp.fwd_kmer_positions + pp.rev_kmer_positions;
            if (stage1_primer_score >= 2.0 &&
                static_cast<uint32_t>(stage1_primer_score) > total_pos) {
                std::fprintf(stderr, "Error: -stage1_primer_score exceeds positions for %s\n",
                             pp.query_id.c_str());
                return 1;
            }
            uint32_t s2 = std::max(pp.fwd_kmer_positions, pp.rev_kmer_positions)
                          + static_cast<uint32_t>(stage2_primer_score_add);
            if (s2 > total_pos) {
                std::fprintf(stderr, "Error: stage2 threshold exceeds positions for %s\n",
                             pp.query_id.c_str());
                return 1;
            }
            min_stage2_score = std::min(min_stage2_score, s2);
        }
        if (stage1_primer_score > 0 && stage1_primer_score <= 1.0) {
            base_req.stage1_min_score_frac_x10000 =
                static_cast<uint16_t>(stage1_primer_score * 10000.0);
        } else if (stage1_primer_score >= 2.0) {
            base_req.stage1_min_score = static_cast<uint16_t>(stage1_primer_score);
        }
        base_req.stage2_min_score = static_cast<uint16_t>(min_stage2_score);
        base_req.has_stage2_min_score = 1;
        if (!cli.has("-stage2_max_gap")) {
            base_req.stage2_max_gap = static_cast<uint16_t>(insert_length);
        }
    }

#ifdef IKAFSSN_ENABLE_HTTP
    if (has_http) {
        install_sigint_handler();
        std::string http_url = cli.get_string("-http");
        std::string err;
        if (!ensure_jobs_root(err)) {
            std::fprintf(stderr, "%s\n", err.c_str());
            return 1;
        }
        std::string root = default_jobs_root();

        // Determine batch size from server info.
        int batch_size = static_cast<int>(queries.size());
        if (server_info.max_nseq_per_req > 0)
            batch_size = std::min(batch_size, server_info.max_nseq_per_req);
        if (batch_size < 1) batch_size = 1;

        GroupMeta gm;
        gm.group_id = make_uuidv4();
        gm.submitted_at =
            std::chrono::duration_cast<std::chrono::seconds>(
                std::chrono::system_clock::now().time_since_epoch()).count();
        gm.httpd_url = http_url;
        gm.db = base_req.db;
        gm.query_file_path_abs = (input_path == "-")
            ? (root + "/" + gm.group_id + "/stdin.fasta")
            : input_path;
        if (input_path == "-") gm.query_file_sha256 = sha256_string(stdin_content);
        else                    gm.query_file_sha256 = sha256_file(input_path);
        gm.max_nseq_per_req = batch_size;
        gm.k             = resolved_k;
        gm.mode          = base_req.mode;
        gm.t             = base_req.t;
        gm.template_type = base_req.template_type;
        gm.output_format = (outfmt == OutputFormat::kJson) ? "json"
                  : (outfmt == OutputFormat::kSam)  ? "sam"
                  : (outfmt == OutputFormat::kBam)  ? "bam"
                  : "tsv";
        gm.output_path = output_path;
        gm.compression_level = compression_level;

        // Persist stdin content as <group>/stdin.fasta so a later -resume can
        // recover the original input even if the user closed their pipe.
        if (input_path == "-") {
            std::error_code ec;
            std::filesystem::create_directories(
                std::filesystem::path(gm.query_file_path_abs).parent_path(), ec);
            std::ofstream sf(gm.query_file_path_abs, std::ios::binary);
            sf.write(stdin_content.data(),
                     static_cast<std::streamsize>(stdin_content.size()));
        }

        for (size_t off = 0; off < queries.size(); off += batch_size) {
            size_t end = std::min(off + static_cast<size_t>(batch_size),
                                   queries.size());
            std::string job_id = make_uuidv4();
            gm.job_ids.push_back(job_id);

            JobMeta jm;
            jm.job_id = job_id;
            jm.group_id = gm.group_id;
            jm.n_seqs = static_cast<int32_t>(end - off);
            jm.fasta_file = gm.query_file_path_abs;
            jm.seq_index_range = {static_cast<int32_t>(off),
                                   static_cast<int32_t>(end)};
            jm.status = "submitted";
            jm.submitted_at = gm.submitted_at;
            std::string werr;
            write_job_meta(root, jm, werr);

            std::vector<std::string> deflines;
            for (size_t i = off; i < end; i++) deflines.push_back(queries[i].id);
            std::string derr;
            write_job_deflines(root, gm.group_id, job_id, deflines, derr);
        }
        write_group_meta(root, gm, err);

        // Submit jobs (best-effort; failures leave the meta on disk so
        // -resume can pick them up).
        for (size_t i = 0; i < gm.job_ids.size(); i++) {
            const auto& job_id = gm.job_ids[i];
            size_t off = i * static_cast<size_t>(batch_size);
            size_t end = std::min(off + static_cast<size_t>(batch_size),
                                   queries.size());
            SearchRequest req = base_req;
            for (size_t j = off; j < end; j++) {
                req.queries.push_back({queries[j].id, queries[j].sequence});
            }
            std::string serr;
            auto rc = http_submit_job(http_url, job_id, req, serr, auth);
            JobMeta jm;
            std::string e;
            read_job_meta(root, gm.group_id, job_id, jm, e);
            if (rc == AsyncHttpOutcome::kOk || rc == AsyncHttpOutcome::kConflict) {
                jm.status = "queued";
            } else {
                jm.status = "submitted"; // will be retried by -resume
                jm.error_message = serr;
                logger.error("submit %s failed: %s", job_id.c_str(), serr.c_str());
            }
            write_job_meta(root, jm, e);
        }

        if (cli.has("-submit_only")) {
            std::printf("%s\n", gm.group_id.c_str());
            return 0;
        }

        PollLoop loop(root, gm.group_id, auth, logger);
        if (!loop.run()) {
            if (PollLoop::sigint_flag().load()) return 130;
        }
        // Re-read group meta to pick up worker-updated counters.
        GroupMeta gm2;
        std::string g2err;
        if (!read_group_meta(root, gm.group_id, gm2, g2err)) {
            logger.error("post-poll group_meta read: %s", g2err.c_str());
            gm2 = gm;
        }
        return finalize_group(root, gm2, logger);
    }
#endif

    // -- socket / TCP path: synchronous loop with checkpoint --
    std::string input_sha256;
    if (input_path == "-") input_sha256 = sha256_string(stdin_content);
    else                    input_sha256 = sha256_file(input_path);
    std::string seqidlist_sha256;
    std::string neg_seqidlist_sha256;
    if (!seqidlist_path.empty())
        seqidlist_sha256 = sha256_file(seqidlist_path);
    if (!neg_seqidlist_path.empty())
        neg_seqidlist_sha256 = sha256_file(neg_seqidlist_path);

    DbStats db_stats = resolve_db_stats(server_info, base_req.db, resolved_k);
    std::string options_text = build_options_text(
        base_req, db_stats, resolved_k, outfmt,
        seqidlist_sha256, neg_seqidlist_sha256);
    if (has_primer) {
        options_text += "primer=1\n";
        options_text += "insert_length=" + std::to_string(insert_length) + "\n";
        options_text += "stage1_primer_score=" + std::to_string(stage1_primer_score) + "\n";
        options_text += "stage2_primer_score_add=" + std::to_string(stage2_primer_score_add) + "\n";
    }

    Checkpoint::Config ckpt_cfg;
    ckpt_cfg.output_path = output_path;
    ckpt_cfg.input_path = input_path;
    ckpt_cfg.ix_name = ix_name;
    ckpt_cfg.resolved_k = resolved_k;
    ckpt_cfg.output_format = outfmt;
    ckpt_cfg.compression_level = compression_level;

    Checkpoint ckpt(ckpt_cfg, logger);
    LockGuard lock;
    ckpt.acquire_lock(lock);

    std::unordered_set<std::string> completed_seqids;
    int next_batch_num = 0;
    if (ckpt.exists()) {
        if (!lock.locked()) ckpt.acquire_lock(lock);
        if (!ckpt.resume(options_text, input_sha256,
                          completed_seqids, next_batch_num)) {
            lock.release();
            ckpt.cleanup();
            if (!ckpt.initialize(options_text, input_sha256, stdin_content)) return 1;
            ckpt.acquire_lock(lock);
        }
    } else {
        if (!ckpt.initialize(options_text, input_sha256, stdin_content)) return 1;
        if (!lock.locked()) ckpt.acquire_lock(lock);
    }

    std::vector<FastaRecord> remaining;
    for (const auto& q : queries) {
        if (completed_seqids.find(q.id) == completed_seqids.end()) {
            remaining.push_back(q);
        }
    }
    std::unordered_map<std::string, std::string> query_map;
    for (const auto& q : queries) query_map[q.id] = q.sequence;

    int batch_size = static_cast<int>(remaining.size());
    if (batch_size == 0) batch_size = 1;
    if (server_info.max_nseq_per_req > 0)
        batch_size = std::min(batch_size, static_cast<int>(server_info.max_nseq_per_req));
    if (server_info.max_queue_size > 0) {
        int available = server_info.max_queue_size - server_info.queue_depth;
        if (available > 0) batch_size = std::min(batch_size, available);
    }
    if (batch_size <= 0) batch_size = 1;

    bool has_skipped = false;
    uint8_t resp_mode = 0;
    uint8_t resp_stage1_score = 0;
    bool resp_stage3_traceback = false;
    // first_response stays true until a response is seen; if there is no
    // remaining work and the checkpoint already has response metadata, the
    // first response was already recorded.
    bool first_response = true;
    if (remaining.empty() && ckpt.read_response_meta(resp_mode, resp_stage1_score,
                                                     resp_stage3_traceback)) {
        first_response = false;
    }

    static constexpr int retry_delays[] = {30, 60, 120};
    static constexpr int num_retry_delays =
        static_cast<int>(sizeof(retry_delays) / sizeof(retry_delays[0]));
    int batch_num = next_batch_num;
    size_t sent = 0;
    while (sent < remaining.size()) {
        size_t batch_end = std::min(sent + static_cast<size_t>(batch_size),
                                     remaining.size());
        std::vector<std::string> batch_seqids;
        for (size_t i = sent; i < batch_end; i++) {
            batch_seqids.push_back(remaining[i].id);
        }
        ckpt.write_batch_seqids(batch_num, batch_seqids);
        SearchRequest req = base_req;
        for (size_t i = sent; i < batch_end; i++) {
            req.queries.push_back({remaining[i].id, remaining[i].sequence});
        }
        sent = batch_end;
        for (int attempt = 0; ; attempt++) {
            SearchResponse resp;
            if (!execute_search(cli, req, resp, logger)) {
                lock.release();
                return 1;
            }
            if (resp.status != 0) {
                std::fprintf(stderr, "Error: server returned status %d\n", resp.status);
                lock.release();
                return 1;
            }
            if (first_response) {
                resp_mode = resp.mode;
                resp_stage1_score = resp.stage1_score;
                resp_stage3_traceback = (resp.stage3_traceback != 0);
                first_response = false;
                ckpt.write_response_meta(resp_mode, resp_stage1_score, resp_stage3_traceback);
            }
            std::vector<OutputHit> batch_hits;
            collect_results(resp, batch_hits, has_skipped);
            if (resp.rejected_qseqids.empty()) {
                ckpt.write_batch_results(batch_num, batch_hits,
                                          resp_mode, resp_stage1_score,
                                          resp_stage3_traceback);
                batch_num++;
                break;
            }
            std::unordered_set<std::string> rejected_set(
                resp.rejected_qseqids.begin(), resp.rejected_qseqids.end());
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
            int delay = retry_delays[std::min(attempt, num_retry_delays - 1)];
            std::this_thread::sleep_for(std::chrono::seconds(delay));
            batch_seqids.clear();
            req = base_req;
            for (const auto& qid : resp.rejected_qseqids) {
                auto it = query_map.find(qid);
                if (it != query_map.end()) {
                    req.queries.push_back({qid, it->second});
                    batch_seqids.push_back(qid);
                }
            }
            if (req.queries.empty()) break;
            ckpt.write_batch_seqids(batch_num, batch_seqids);
        }
    }
    if (remaining.empty() && first_response) {
        if (!ckpt.read_response_meta(resp_mode, resp_stage1_score,
                                      resp_stage3_traceback)) {
            lock.release();
            return 1;
        }
    }
    if (!ckpt.merge_results(output_path, resp_mode, resp_stage1_score,
                             resp_stage3_traceback)) {
        lock.release();
        return 1;
    }
    lock.release();
    ckpt.cleanup();
    return has_skipped ? 2 : 0;
}
