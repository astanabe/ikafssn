#include "ikafssnserver/server.hpp"
#include "core/spaced_seed.hpp"
#include "core/version.hpp"
#include "util/cli_parser.hpp"
#include "util/cli_validators.hpp"
#include "util/common_init.hpp"
#include "util/simd_dispatch.hpp"
#include "util/context_parser.hpp"
#include "util/size_parser.hpp"
#include "io/volume_discovery.hpp"

#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <set>

using namespace ikafssn;

static Server* g_server = nullptr;

static void signal_handler(int /*sig*/) {
    if (g_server) {
        g_server->request_shutdown();
    }
}

static void print_usage(const char* prog) {
    print_version_header("ikafssnserver");
    std::fprintf(stderr,
        "Usage: %s [options]\n"
        "\n"
        "Required:\n"
        "  -ix <prefix>             Index prefix (repeatable for multi-DB)\n"
        "\n"
        "Listener (at least one required):\n"
        "  -socket <path>           UNIX domain socket path\n"
        "  -tcp <host>:<port>       TCP listen address\n"
        "\n"
        "Index-variant load filter (each unset field is a wildcard):\n"
        "  -k <int>                 Load only this k-mer length (default: any)\n"
        "  -t <int>                 Load only this template length: 0=contiguous, 13/15/16/18/21\n"
        "                           (default: any; -t 0 loads the contiguous index only)\n"
        "  -template_type <string>  coding, optimal, contiguous, both (default: any; invalid with -t 0)\n"
        "  -min_seq_length <int>    Load only variants with this min_seq_length (default: any)\n"
        "  -min_length_split <int>  Load only variants with this min_length_split (default: any)\n"
        "  -overlap_length <int>    Load only variants with this overlap_length (default: any)\n"
        "  -max_freq_build <int>    Load only variants with this max_freq_build (default: any)\n"
        "  -max_degen_expand <int>  Load only variants with this max_degen_expand (default: any).\n"
        "                           This is a load filter; the query-side max_degen_expand is\n"
        "                           supplied by the client (server fallback when omitted: 16).\n"
        "\n"
        "Options:\n"
        "  -nthread <int>           Worker threads (default: all cores)\n"
        "  -max_queue_size <int>    Max concurrent query sequences globally (default: 1024)\n"
        "  -max_nseq_per_req <int>  Max sequences accepted per request (default: thread count)\n"
        "  -max_concurrent_search <int>  Cap concurrent searches sharing -memory_limit's\n"
        "                           residual posting_budget. 0 = unlimited (default).\n"
        "                           When N >= 1, at most N searches run at once, bounding\n"
        "                           in-flight posting heap (Stage 3) to posting_budget.\n"
        "  -pid <path>              PID file path\n"
        "  -db <path>               BLAST DB path for mode 3 (repeatable, paired with -ix;\n"
        "                           default: same as corresponding -ix prefix)\n"
        "  -stage2_min_score <int>  Default minimum chain score (default: 0 = adaptive)\n"
        "  -stage2_max_gap <int>    Default chaining gap tolerance (default: 100)\n"
        "  -stage2_max_lookback <int>  Default chaining DP lookback window (default: 64, 0=unlimited)\n"
        "  -stage2_max_nhit_per_subject <int>  Default max chains per subject (default: 1, 0=unlimited)\n"
        "  -stage2_min_nhit_diag <int>  Default diagonal filter min hits (default: 1)\n"
        "  -stage1_topn <int>       Default Stage 1 candidate limit (default: 0)\n"
        "  -stage1_min_score <num>  Default Stage 1 minimum score; integer or 0<P<1 fraction (default: 0.5)\n"
        "  -nresult <int>           Default max results per query (default: 0)\n"
        "  -accept_qdegen <0|1>     Default accept queries with degenerate bases (default: 1)\n"
        "  -context_extend <value>  Default context extension (int=bases, decimal=ratio, default: 2.0)\n"
        "  -stage3_traceback <0|1>  Default traceback mode (default: 0)\n"
        "  -stage3_gapopen <int>    Default gap open penalty (default: 10)\n"
        "  -stage3_gapext <int>     Default gap extension penalty (default: 1)\n"
        "  -stage3_min_ppositive <num> Default min percent positive (default: 0)\n"
        "  -stage3_min_npositive <int> Default min positive-scoring positions (default: 0)\n"
        "  -stage3_score_matrix <name> Default score matrix: degmatch, dnafull, nuc44 (default: degmatch)\n"
        "  -memory_limit <size>     Search memory budget (default: half of RAM)\n"
        "                           Pins .khx/.ksx metadata; residual caps the\n"
        "                           Stage 3 posting heap and concurrent-search pool.\n"
        "                           Accepts K, M, G suffixes\n"
        "  -shutdown_timeout <int>  Graceful shutdown timeout in seconds (default: 30)\n"
        "  -v, --verbose            Verbose logging\n",
        prog);
}

int main(int argc, char* argv[]) {
    CliParser cli(argc, argv);

    if (check_version(cli, "ikafssnserver")) return 0;

    if (cli.has("-h") || cli.has("-help")) {
        print_usage(argv[0]);
        return 0;
    }

    if (!cli.has("-ix")) {
        std::fprintf(stderr, "Error: -ix is required\n");
        print_usage(argv[0]);
        return 1;
    }

    if (!cli.has("-socket") && !cli.has("-tcp")) {
        std::fprintf(stderr, "Error: at least one of -socket or -tcp is required\n");
        print_usage(argv[0]);
        return 1;
    }

    // -t 0 with an explicit -template_type is contradictory.
    {
        std::string err;
        if (!validate_t_template_type_combo(cli, err)) {
            std::fprintf(stderr, "%s\n", err.c_str());
            return 1;
        }
    }

    // Build the index-variant load filter (ikafssninfo-local semantics).
    VariantFilter load_filter;
    {
        load_filter.has_k = cli.has("-k");
        load_filter.k = cli.get_int("-k", 0);
        load_filter.t_is_wildcard = !cli.has("-t");
        if (cli.has("-t")) {
            std::string err;
            if (!parse_spaced_seed_t(cli, load_filter.t, err)) {
                std::fprintf(stderr, "%s\n", err.c_str());
                return 1;
            }
        }
        if (cli.has("-template_type")) {
            const std::string s = cli.get_string("-template_type");
            if (s != "coding" && s != "optimal" && s != "contiguous" &&
                s != "both" && s != "coding_and_optimal") {
                std::fprintf(stderr, "Error: invalid -template_type '%s' "
                             "(coding, optimal, contiguous, both)\n", s.c_str());
                return 1;
            }
            load_filter.has_template_type = true;
            load_filter.template_type =
                static_cast<uint8_t>(template_type_from_string(s));
        }
        load_filter.has_min_seq_length = cli.has("-min_seq_length");
        load_filter.min_seq_length =
            static_cast<uint32_t>(cli.get_int("-min_seq_length", 0));
        load_filter.has_min_length_split = cli.has("-min_length_split");
        load_filter.min_length_split =
            static_cast<uint32_t>(cli.get_int("-min_length_split", 0));
        load_filter.has_overlap_length = cli.has("-overlap_length");
        load_filter.overlap_length =
            static_cast<uint32_t>(cli.get_int("-overlap_length", 0));
        load_filter.has_max_freq_build = cli.has("-max_freq_build");
        load_filter.max_freq_build = cli.has("-max_freq_build")
            ? std::strtoull(cli.get_string("-max_freq_build").c_str(), nullptr, 10)
            : 1;
        load_filter.has_max_degen_expand = cli.has("-max_degen_expand");
        if (cli.has("-max_degen_expand")) {
            int v = cli.get_int("-max_degen_expand", 0);
            if (v < 0 || v > 256) {
                std::fprintf(stderr,
                    "Error: -max_degen_expand must be between 0 and 256\n");
                return 1;
            }
            load_filter.max_degen_expand = static_cast<uint32_t>(v);
        }
    }

    // Build db_entries from -ix and -db lists
    auto ix_list = cli.get_strings("-ix");
    auto db_list = cli.get_strings("-db");

    if (!db_list.empty() && db_list.size() != ix_list.size()) {
        std::fprintf(stderr,
            "Error: -db count (%zu) must be 0 or equal to -ix count (%zu)\n",
            db_list.size(), ix_list.size());
        return 1;
    }

    // Check for duplicate DB names
    {
        std::set<std::string> seen_names;
        for (const auto& ix : ix_list) {
            auto parts = parse_index_prefix(ix);
            if (!seen_names.insert(parts.db).second) {
                std::fprintf(stderr,
                    "Error: duplicate database name '%s' (from -ix %s)\n",
                    parts.db.c_str(), ix.c_str());
                return 1;
            }
        }
    }

    ServerConfig config;
    config.variant_filter = load_filter;
    for (size_t i = 0; i < ix_list.size(); i++) {
        ServerConfig::DbEntry entry;
        entry.ix_prefix = ix_list[i];
        entry.db_path = (i < db_list.size()) ? db_list[i] : ix_list[i];
        config.db_entries.push_back(std::move(entry));
    }

    config.unix_socket_path = cli.get_string("-socket");
    config.tcp_addr = cli.get_string("-tcp");
    config.pid_file = cli.get_string("-pid");
    config.nthread = cli.get_int("-nthread", 0);
    config.max_queue_size = cli.get_int("-max_queue_size", 0);
    if (config.max_queue_size < 0) {
        std::fprintf(stderr, "Error: -max_queue_size must be >= 0\n");
        return 1;
    }
    config.max_nseq_per_req = cli.get_int("-max_nseq_per_req", 0);
    if (config.max_nseq_per_req < 0) {
        std::fprintf(stderr, "Error: -max_nseq_per_req must be >= 0\n");
        return 1;
    }
    config.max_concurrent_search = cli.get_int("-max_concurrent_search", 0);
    if (config.max_concurrent_search < 0) {
        std::fprintf(stderr, "Error: -max_concurrent_search must be >= 0\n");
        return 1;
    }
    config.shutdown_timeout = cli.get_int("-shutdown_timeout", 30);

    if (cli.has("-memory_limit")) {
        std::string mem_str = cli.get_string("-memory_limit");
        config.memory_limit = parse_size_string(mem_str);
        if (config.memory_limit == 0) {
            std::fprintf(stderr, "Error: invalid -memory_limit '%s'\n", mem_str.c_str());
            return 1;
        }
    }

    Logger logger = make_logger(cli);
    init_simd_dispatch(&logger);
    config.log_level = logger.level();

    // Search config
    config.search_config.stage1.stage1_topn =
        static_cast<uint32_t>(cli.get_int("-stage1_topn", 0));
    {
        double min_s1 = cli.get_double("-stage1_min_score", 0.5);
        if (min_s1 > 0 && min_s1 < 1.0) {
            config.search_config.min_stage1_score_frac = min_s1;
        } else {
            config.search_config.stage1.min_stage1_score = static_cast<uint32_t>(min_s1);
        }
    }
    config.search_config.stage2.max_gap =
        static_cast<uint32_t>(cli.get_int("-stage2_max_gap", 100));
    config.search_config.stage2.chain_max_lookback =
        static_cast<uint32_t>(cli.get_int("-stage2_max_lookback", 64));
    config.search_config.stage2.max_nhit_per_subject =
        static_cast<uint32_t>(cli.get_int("-stage2_max_nhit_per_subject", 1));
    config.search_config.stage2.min_nhit_diag =
        static_cast<uint32_t>(cli.get_int("-stage2_min_nhit_diag", 1));
    config.search_config.stage2.min_score =
        static_cast<uint32_t>(cli.get_int("-stage2_min_score", 0));
    config.search_config.nresult =
        static_cast<uint32_t>(cli.get_int("-nresult", 0));
    config.search_config.accept_qdegen =
        static_cast<uint8_t>(cli.get_int("-accept_qdegen", 1));
    // Query-side max_degen_expand is a client-supplied parameter; the server
    // only needs an internal fallback (16, matching the client default) for
    // requests that omit it.  -max_degen_expand is an index load filter
    // (handled above), not the query default.
    config.search_config.max_degen_expand = 16;

    // Stage 3 config
    config.stage3_config.gapopen = cli.get_int("-stage3_gapopen", 10);
    config.stage3_config.gapext = cli.get_int("-stage3_gapext", 1);
    config.stage3_config.traceback = (cli.get_int("-stage3_traceback", 0) != 0);
    config.stage3_config.min_ppositive = cli.get_double("-stage3_min_ppositive", 0.0);
    config.stage3_config.min_npositive = static_cast<uint32_t>(cli.get_int("-stage3_min_npositive", 0));
    {
        std::string err;
        if (!parse_score_matrix(cli, config.stage3_config.score_matrix, err)) {
            std::fprintf(stderr, "%s\n", err.c_str());
            return 1;
        }
    }

    // Parse -context_extend: integer (bases) or decimal (query length multiplier)
    {
        ContextParam ctx_param;
        std::string err;
        if (!parse_context(cli.get_string("-context_extend", "2.0"), ctx_param, err)) {
            std::fprintf(stderr, "%s\n", err.c_str());
            return 1;
        }
        config.context_is_ratio = ctx_param.is_ratio;
        config.context_ratio = ctx_param.ratio;
        config.context_abs = ctx_param.abs;
    }

    Server server;
    g_server = &server;

    // Install signal handlers
    struct sigaction sa;
    sa.sa_handler = signal_handler;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = 0;
    sigaction(SIGTERM, &sa, nullptr);
    sigaction(SIGINT, &sa, nullptr);

    int ret = server.run(config);
    g_server = nullptr;
    return ret;
}
