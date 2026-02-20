#include "ikafssnclient/socket_client.hpp"
#include "io/fasta_reader.hpp"
#include "io/seqidlist_reader.hpp"
#include "io/result_writer.hpp"
#include "protocol/messages.hpp"
#include "util/cli_parser.hpp"
#include "util/socket_utils.hpp"
#include "util/logger.hpp"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

using namespace ikafssn;

static void print_usage(const char* prog) {
    std::fprintf(stderr,
        "Usage: %s [options]\n"
        "\n"
        "Connection (one required):\n"
        "  -socket <path>           UNIX domain socket path\n"
        "  -tcp <host>:<port>       TCP server address\n"
        "\n"
        "Required:\n"
        "  -query <path>            Query FASTA file (- for stdin)\n"
        "\n"
        "Options:\n"
        "  -o <path>                Output file (default: stdout)\n"
        "  -k <int>                 K-mer size (default: server default)\n"
        "  -min_score <int>         Minimum chain score (default: server default)\n"
        "  -max_gap <int>           Chaining gap tolerance (default: server default)\n"
        "  -max_freq <int>          High-freq k-mer skip threshold (default: server default)\n"
        "  -min_diag_hits <int>     Diagonal filter min hits (default: server default)\n"
        "  -stage1_topn <int>       Stage 1 candidate limit (default: server default)\n"
        "  -min_stage1_score <int>  Stage 1 minimum score (default: server default)\n"
        "  -num_results <int>       Max results per query (default: server default)\n"
        "  -seqidlist <path>        Include only listed accessions\n"
        "  -negative_seqidlist <path>  Exclude listed accessions\n"
        "  -outfmt <tab|json>       Output format (default: tab)\n"
        "  -v, --verbose            Verbose logging\n",
        prog);
}

int main(int argc, char* argv[]) {
    CliParser cli(argc, argv);

    if (cli.has("-h") || cli.has("--help")) {
        print_usage(argv[0]);
        return 0;
    }

    if (!cli.has("-socket") && !cli.has("-tcp")) {
        std::fprintf(stderr, "Error: one of -socket or -tcp is required\n");
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

    // Build search request
    SearchRequest req;
    req.k = static_cast<uint8_t>(cli.get_int("-k", 0));
    req.min_score = static_cast<uint16_t>(cli.get_int("-min_score", 0));
    req.max_gap = static_cast<uint16_t>(cli.get_int("-max_gap", 0));
    req.max_freq = static_cast<uint32_t>(cli.get_int("-max_freq", 0));
    req.min_diag_hits = static_cast<uint8_t>(cli.get_int("-min_diag_hits", 0));
    req.stage1_topn = static_cast<uint16_t>(cli.get_int("-stage1_topn", 0));
    req.min_stage1_score = static_cast<uint16_t>(cli.get_int("-min_stage1_score", 0));
    req.num_results = static_cast<uint16_t>(cli.get_int("-num_results", 0));

    // Seqidlist
    if (cli.has("-seqidlist")) {
        req.seqidlist_mode = SeqidlistMode::kInclude;
        req.seqids = read_seqidlist(cli.get_string("-seqidlist"));
        logger.info("Loaded %zu accessions (include mode)", req.seqids.size());
    } else if (cli.has("-negative_seqidlist")) {
        req.seqidlist_mode = SeqidlistMode::kExclude;
        req.seqids = read_seqidlist(cli.get_string("-negative_seqidlist"));
        logger.info("Loaded %zu accessions (exclude mode)", req.seqids.size());
    }

    // Add queries
    for (const auto& q : queries) {
        req.queries.push_back({q.id, q.sequence});
    }

    // Connect to server
    int fd = -1;
    if (cli.has("-socket")) {
        std::string sock_path = cli.get_string("-socket");
        fd = unix_connect(sock_path);
        if (fd < 0) {
            std::fprintf(stderr, "Error: cannot connect to UNIX socket %s\n", sock_path.c_str());
            return 1;
        }
        logger.debug("Connected to UNIX socket %s", sock_path.c_str());
    } else {
        std::string tcp_addr = cli.get_string("-tcp");
        fd = tcp_connect(tcp_addr);
        if (fd < 0) {
            std::fprintf(stderr, "Error: cannot connect to TCP %s\n", tcp_addr.c_str());
            return 1;
        }
        logger.debug("Connected to TCP %s", tcp_addr.c_str());
    }

    // Send request and receive response
    SearchResponse resp;
    if (!socket_search(fd, req, resp)) {
        std::fprintf(stderr, "Error: search request failed\n");
        close_fd(fd);
        return 1;
    }
    close_fd(fd);

    if (resp.status != 0) {
        std::fprintf(stderr, "Error: server returned status %d\n", resp.status);
        return 1;
    }

    logger.info("Received response: k=%d, %zu query result(s)", resp.k, resp.results.size());

    // Convert response to OutputHit format
    std::vector<OutputHit> all_hits;
    for (const auto& qr : resp.results) {
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
            oh.volume = hit.volume;
            all_hits.push_back(oh);
        }
    }

    // Write output
    if (output_path.empty()) {
        write_results(std::cout, all_hits, outfmt);
    } else {
        std::ofstream out(output_path);
        if (!out.is_open()) {
            std::fprintf(stderr, "Error: cannot open output file %s\n", output_path.c_str());
            return 1;
        }
        write_results(out, all_hits, outfmt);
    }

    logger.info("Done. %zu hit(s) reported.", all_hits.size());
    return 0;
}
