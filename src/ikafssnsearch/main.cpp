#include "core/config.hpp"
#include "core/types.hpp"
#include "core/kmer_encoding.hpp"
#include "index/kix_reader.hpp"
#include "index/kpx_reader.hpp"
#include "index/ksx_reader.hpp"
#include "search/oid_filter.hpp"
#include "search/volume_searcher.hpp"
#include "io/fasta_reader.hpp"
#include "io/seqidlist_reader.hpp"
#include "io/result_writer.hpp"
#include "util/cli_parser.hpp"
#include "util/logger.hpp"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <regex>
#include <string>
#include <vector>

using namespace ikafssn;

static void print_usage(const char* prog) {
    std::fprintf(stderr,
        "Usage: %s [options]\n"
        "\n"
        "Required:\n"
        "  -ix <dir>                Index directory\n"
        "  -query <path>            Query FASTA file (- for stdin)\n"
        "\n"
        "Options:\n"
        "  -o <path>                Output file (default: stdout)\n"
        "  -min_score <int>         Minimum chain score (default: 3)\n"
        "  -max_gap <int>           Chaining diagonal gap tolerance (default: 100)\n"
        "  -max_freq <int>          High-frequency k-mer skip threshold (default: auto)\n"
        "  -min_diag_hits <int>     Diagonal filter min hits (default: 2)\n"
        "  -stage1_topn <int>       Stage 1 candidate limit (default: 500)\n"
        "  -min_stage1_score <int>  Stage 1 minimum score (default: 2)\n"
        "  -num_results <int>       Max results per query (default: 50)\n"
        "  -seqidlist <path>        Include only listed accessions\n"
        "  -negative_seqidlist <path>  Exclude listed accessions\n"
        "  -outfmt <tab|json>       Output format (default: tab)\n"
        "  -v, --verbose            Verbose logging\n",
        prog);
}

struct VolumeFiles {
    std::string kix_path;
    std::string kpx_path;
    std::string ksx_path;
    uint16_t volume_index;
    int k;
};

static std::vector<VolumeFiles> discover_volumes(const std::string& ix_dir) {
    std::vector<VolumeFiles> volumes;
    std::regex kix_pattern(R"((.+)\.(\d+)\.(\d+)mer\.kix)");

    for (const auto& entry : std::filesystem::directory_iterator(ix_dir)) {
        if (!entry.is_regular_file()) continue;
        std::string fname = entry.path().filename().string();
        std::smatch m;
        if (std::regex_match(fname, m, kix_pattern)) {
            VolumeFiles vf;
            vf.volume_index = static_cast<uint16_t>(std::stoi(m[2].str()));
            vf.k = std::stoi(m[3].str());
            std::string base = ix_dir + "/" + m[1].str() + "." + m[2].str() + "." + m[3].str() + "mer";
            vf.kix_path = base + ".kix";
            vf.kpx_path = base + ".kpx";
            vf.ksx_path = base + ".ksx";
            volumes.push_back(vf);
        }
    }

    std::sort(volumes.begin(), volumes.end(),
              [](const VolumeFiles& a, const VolumeFiles& b) {
                  return a.volume_index < b.volume_index;
              });
    return volumes;
}

int main(int argc, char* argv[]) {
    CliParser cli(argc, argv);

    if (cli.has("-h") || cli.has("--help")) {
        print_usage(argv[0]);
        return 0;
    }

    if (!cli.has("-ix") || !cli.has("-query")) {
        print_usage(argv[0]);
        return 1;
    }

    // Check mutually exclusive options
    if (cli.has("-seqidlist") && cli.has("-negative_seqidlist")) {
        std::fprintf(stderr, "Error: -seqidlist and -negative_seqidlist are mutually exclusive\n");
        return 1;
    }

    std::string ix_dir = cli.get_string("-ix");
    std::string query_path = cli.get_string("-query");
    std::string output_path = cli.get_string("-o");
    bool verbose = cli.has("-v") || cli.has("--verbose");

    Logger logger(verbose ? Logger::kDebug : Logger::kInfo);

    // Search config
    SearchConfig config;
    config.stage1.max_freq = static_cast<uint32_t>(cli.get_int("-max_freq", 0));
    config.stage1.stage1_topn = static_cast<uint32_t>(cli.get_int("-stage1_topn", 500));
    config.stage1.min_stage1_score = static_cast<uint32_t>(cli.get_int("-min_stage1_score", 2));
    config.stage2.max_gap = static_cast<uint32_t>(cli.get_int("-max_gap", 100));
    config.stage2.min_diag_hits = static_cast<uint32_t>(cli.get_int("-min_diag_hits", 2));
    config.stage2.min_score = static_cast<uint32_t>(cli.get_int("-min_score", 3));
    config.num_results = static_cast<uint32_t>(cli.get_int("-num_results", 50));

    // Output format
    OutputFormat outfmt = OutputFormat::kTab;
    std::string outfmt_str = cli.get_string("-outfmt", "tab");
    if (outfmt_str == "json") {
        outfmt = OutputFormat::kJson;
    } else if (outfmt_str != "tab") {
        std::fprintf(stderr, "Error: unknown output format '%s'\n", outfmt_str.c_str());
        return 1;
    }

    // Discover volumes
    auto volumes = discover_volumes(ix_dir);
    if (volumes.empty()) {
        std::fprintf(stderr, "Error: no index files found in %s\n", ix_dir.c_str());
        return 1;
    }

    int k = volumes[0].k;
    logger.info("Found %zu volume(s), k=%d", volumes.size(), k);

    // Read query FASTA
    auto queries = read_fasta(query_path);
    if (queries.empty()) {
        std::fprintf(stderr, "Error: no query sequences found\n");
        return 1;
    }
    logger.info("Read %zu query sequence(s)", queries.size());

    // Read seqidlist if specified
    std::vector<std::string> seqidlist;
    OidFilterMode filter_mode = OidFilterMode::kNone;
    if (cli.has("-seqidlist")) {
        seqidlist = read_seqidlist(cli.get_string("-seqidlist"));
        filter_mode = OidFilterMode::kInclude;
        logger.info("Loaded %zu accessions from seqidlist (include mode)", seqidlist.size());
    } else if (cli.has("-negative_seqidlist")) {
        seqidlist = read_seqidlist(cli.get_string("-negative_seqidlist"));
        filter_mode = OidFilterMode::kExclude;
        logger.info("Loaded %zu accessions from seqidlist (exclude mode)", seqidlist.size());
    }

    // Collect all output hits
    std::vector<OutputHit> all_hits;

    // Process each volume
    for (const auto& vol : volumes) {
        logger.info("Searching volume %u...", vol.volume_index);

        KixReader kix;
        if (!kix.open(vol.kix_path)) {
            std::fprintf(stderr, "Error: cannot open %s\n", vol.kix_path.c_str());
            return 1;
        }
        KpxReader kpx;
        if (!kpx.open(vol.kpx_path)) {
            std::fprintf(stderr, "Error: cannot open %s\n", vol.kpx_path.c_str());
            return 1;
        }
        KsxReader ksx;
        if (!ksx.open(vol.ksx_path)) {
            std::fprintf(stderr, "Error: cannot open %s\n", vol.ksx_path.c_str());
            return 1;
        }

        // Build OID filter for this volume
        OidFilter filter;
        if (filter_mode != OidFilterMode::kNone) {
            filter.build(seqidlist, ksx, filter_mode);
        }

        // Search each query
        for (const auto& query : queries) {
            SearchResult sr;
            if (k < K_TYPE_THRESHOLD) {
                sr = search_volume<uint16_t>(
                    query.id, query.sequence, k, kix, kpx, ksx, filter, config);
            } else {
                sr = search_volume<uint32_t>(
                    query.id, query.sequence, k, kix, kpx, ksx, filter, config);
            }

            // Convert to OutputHit
            for (const auto& cr : sr.hits) {
                OutputHit oh;
                oh.query_id = sr.query_id;
                oh.accession = std::string(ksx.accession(cr.seq_id));
                oh.strand = cr.is_reverse ? '-' : '+';
                oh.q_start = cr.q_start;
                oh.q_end = cr.q_end;
                oh.s_start = cr.s_start;
                oh.s_end = cr.s_end;
                oh.score = cr.score;
                oh.volume = vol.volume_index;
                all_hits.push_back(oh);
            }
        }
    }

    // Sort final results by (query_id, score desc)
    std::sort(all_hits.begin(), all_hits.end(),
              [](const OutputHit& a, const OutputHit& b) {
                  if (a.query_id != b.query_id) return a.query_id < b.query_id;
                  return a.score > b.score;
              });

    // Truncate per query
    {
        std::vector<OutputHit> truncated;
        std::string cur_qid;
        uint32_t cur_count = 0;
        for (const auto& h : all_hits) {
            if (h.query_id != cur_qid) {
                cur_qid = h.query_id;
                cur_count = 0;
            }
            if (cur_count < config.num_results) {
                truncated.push_back(h);
                cur_count++;
            }
        }
        all_hits = std::move(truncated);
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
