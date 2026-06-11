#include "io/compressed_stream.hpp"
#include "io/result_reader.hpp"
#include "io/result_writer.hpp"
#include "ikafssnretrieve/local_retriever.hpp"
#include "core/version.hpp"
#include "util/cli_parser.hpp"
#include "util/common_init.hpp"
#include "util/simd_dispatch.hpp"
#include "util/context_parser.hpp"
#include "util/logger.hpp"

#ifdef IKAFSSN_ENABLE_REMOTE
#include "ikafssnretrieve/efetch_retriever.hpp"
#endif

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

using namespace ikafssn;

static void print_usage(const char* prog) {
    print_version_header("ikafssnretrieve");
    std::fprintf(stderr,
        "Usage: %s [options]\n"
        "\n"
        "Sequence source (one required):\n"
        "  -db <path>              Local BLAST DB prefix\n"
#ifdef IKAFSSN_ENABLE_REMOTE
        "  -remote                 Retrieve from NCBI efetch\n"
#endif
        "\n"
        "Input:\n"
        "  -tsv <path>             Search results file (TSV format)\n"
        "  (none)                  Read from stdin\n"
        "\n"
        "Common options:\n"
        "  -o <path>               Output FASTA file (default: stdout)\n"
        "  -context_extend <value> Context extension: integer=bases, decimal=multiplier of q_len (default: 2.0)\n"
        "  -compression_level <int> Output compression level (default per codec: gzip=6, bzip2=9, xz=6, zstd=3)\n"
        "  -v, --verbose           Verbose logging\n"
#ifdef IKAFSSN_ENABLE_REMOTE
        "\n"
        "Remote options (-remote):\n"
        "  -api_key <key>          NCBI API key (or NCBI_API_KEY env var)\n"
        "  -batch_size <int>       Accessions per batch (default: 100)\n"
        "  -max_nretry <int>       Max retries (default: 3)\n"
        "  -timeout <int>          Request timeout in seconds (default: 30)\n"
        "  -range_threshold <int>  Seq length for individual fetch (default: 100000)\n"
#endif
        ,
        prog);
}

int main(int argc, char* argv[]) {
    CliParser cli(argc, argv);

    if (check_version(cli, "ikafssnretrieve")) return 0;

    if (cli.has("-h") || cli.has("-help")) {
        print_usage(argv[0]);
        return 0;
    }

    bool has_db = cli.has("-db");
    bool has_remote = cli.has("-remote");

    if (!has_db && !has_remote) {
        std::fprintf(stderr, "Error: either -db or -remote is required\n");
        print_usage(argv[0]);
        return 1;
    }
    if (has_db && has_remote) {
        std::fprintf(stderr, "Error: -db and -remote are mutually exclusive\n");
        return 1;
    }

#ifndef IKAFSSN_ENABLE_REMOTE
    if (has_remote) {
        std::fprintf(stderr, "Error: -remote is not available (built without ENABLE_REMOTE_RETRIEVE)\n");
        return 1;
    }
#endif

    Logger logger = make_logger(cli);
    init_simd_dispatch(&logger);

    // Parse -context_extend: integer (bases) or decimal (query length multiplier)
    ContextParam ctx_param;
    {
        std::string err;
        if (!parse_context(cli.get_string("-context_extend", "2.0"), ctx_param, err)) {
            std::fprintf(stderr, "%s\n", err.c_str());
            return 1;
        }
    }

    // Read search results
    std::string tsv_path = cli.get_string("-tsv", "-");
    logger.info("Reading search results from %s",
                tsv_path == "-" ? "stdin" : tsv_path.c_str());

    auto hits = read_results_tsv(tsv_path);
    if (hits.empty()) {
        std::fprintf(stderr, "Error: no valid search results found\n");
        return 1;
    }
    logger.info("Read %zu hit(s)", hits.size());

    // Open output
    std::string output_path = cli.get_string("-o");
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
    std::string open_err;
    auto out_owned = open_output_compressed(output_path, compression_level,
                                              open_err);
    if (!out_owned) {
        std::fprintf(stderr, "%s\n", open_err.c_str());
        return 1;
    }
    std::ostream* out_ptr = out_owned.stream.get();

    uint32_t retrieved = 0;

    if (has_db) {
        // Local retrieval
        std::string db_path = cli.get_string("-db");

        RetrieveOptions opts;
        opts.is_ratio = ctx_param.is_ratio;
        opts.ratio    = ctx_param.ratio;
        opts.context  = ctx_param.abs;
        if (ctx_param.is_ratio) {
            logger.info("Retrieving from local BLAST DB: %s (context ratio=%.4f)",
                        db_path.c_str(), ctx_param.ratio);
        } else {
            logger.info("Retrieving from local BLAST DB: %s", db_path.c_str());
        }
        retrieved = retrieve_local(hits, db_path, opts, *out_ptr);
    }
#ifdef IKAFSSN_ENABLE_REMOTE
    else if (has_remote) {
        // Remote retrieval via NCBI efetch
        EfetchOptions opts;
        opts.is_ratio = ctx_param.is_ratio;
        opts.ratio    = ctx_param.ratio;
        opts.context  = ctx_param.abs;
        if (ctx_param.is_ratio) {
            logger.info("Remote retrieval with context ratio=%.4f (per-hit)",
                        ctx_param.ratio);
        }
        opts.batch_size = static_cast<uint32_t>(cli.get_int("-batch_size", 100));
        opts.max_nretry = static_cast<uint32_t>(cli.get_int("-max_nretry", 3));
        opts.timeout_sec = static_cast<uint32_t>(cli.get_int("-timeout", 30));
        opts.range_threshold = static_cast<uint32_t>(cli.get_int("-range_threshold", 100000));

        // API key: CLI > environment variable
        opts.api_key = cli.get_string("-api_key");
        if (opts.api_key.empty()) {
            const char* env_key = std::getenv("NCBI_API_KEY");
            if (env_key) opts.api_key = env_key;
        }

        logger.info("Retrieving from NCBI efetch (batch_size=%u, range_threshold=%u)",
                     opts.batch_size, opts.range_threshold);
        retrieved = retrieve_remote(hits, opts, *out_ptr);
    }
#endif

    logger.info("Done. %u sequence(s) retrieved.", retrieved);
    return 0;
}
