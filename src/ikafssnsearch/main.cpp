#include "core/config.hpp"
#include "core/types.hpp"
#include "core/version.hpp"
#include "core/kmer_encoding.hpp"
#include "index/kix_reader.hpp"
#include "index/kpx_reader.hpp"
#include "index/ksx_reader.hpp"
#include "index/khx_reader.hpp"
#include "search/oid_filter.hpp"
#include "search/volume_searcher.hpp"
#include "io/fasta_reader.hpp"
#include "io/seqidlist_reader.hpp"
#include "io/result_writer.hpp"
#include "util/cli_parser.hpp"
#include "util/logger.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <mutex>
#include <regex>
#include <set>
#include <string>
#include <thread>
#include <vector>

#include <tbb/parallel_for_each.h>
#include <tbb/task_arena.h>

using namespace ikafssn;

static void print_usage(const char* prog) {
    std::fprintf(stderr,
        "Usage: %s [options]\n"
        "\n"
        "Required:\n"
        "  -ix <prefix>             Index prefix (like blastn -db)\n"
        "  -query <path>            Query FASTA file (- for stdin)\n"
        "\n"
        "Options:\n"
        "  -k <int>                 K-mer size to use (required if multiple k values exist)\n"
        "  -o <path>                Output file (default: stdout)\n"
        "  -threads <int>           Parallel search threads (default: all cores)\n"
        "  -mode <1|2>              1=Stage1 only, 2=Stage1+Stage2 (default: 2)\n"
        "  -stage1_score <1|2>      1=coverscore, 2=matchscore (default: 1)\n"
        "  -sort_score <1|2>        1=stage1 score, 2=chainscore (default: 2)\n"
        "  -min_score <int>         Minimum score (default: 1)\n"
        "  -max_gap <int>           Chaining diagonal gap tolerance (default: 100)\n"
        "  -max_freq <num>          High-frequency k-mer skip threshold (default: 0.5)\n"
        "                           0 < x < 1: fraction of total NSEQ across all volumes\n"
        "                           >= 1: absolute count threshold; 0 = auto\n"
        "  -min_diag_hits <int>     Diagonal filter min hits (default: 2)\n"
        "  -stage1_topn <int>       Stage 1 candidate limit, 0=unlimited (default: 0)\n"
        "  -min_stage1_score <num>  Stage 1 minimum score; integer or 0<P<1 fraction (default: 0.5)\n"
        "  -num_results <int>       Max results per query, 0=unlimited (default: 0)\n"
        "  -seqidlist <path>        Include only listed accessions\n"
        "  -negative_seqidlist <path>  Exclude listed accessions\n"
        "  -strand <-1|1|2>         Strand: 1=plus, -1=minus, 2=both (default: 2)\n"
        "  -accept_qdegen <0|1>     Accept queries with degenerate bases (default: 0)\n"
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

static std::vector<VolumeFiles> discover_volumes(const std::string& ix_prefix, int filter_k = 0) {
    std::vector<VolumeFiles> volumes;

    std::filesystem::path prefix_path(ix_prefix);
    std::string parent_dir = prefix_path.parent_path().string();
    std::string db_name = prefix_path.filename().string();
    if (parent_dir.empty()) parent_dir = ".";

    // Match files: <db_name>.<vol_idx>.<k>mer.kix
    std::regex suffix_pattern(R"((\d+)\.(\d+)mer\.kix)");

    std::string prefix_dot = db_name + ".";
    for (const auto& entry : std::filesystem::directory_iterator(parent_dir)) {
        if (!entry.is_regular_file()) continue;
        std::string fname = entry.path().filename().string();
        if (fname.size() <= prefix_dot.size() || fname.compare(0, prefix_dot.size(), prefix_dot) != 0)
            continue;
        std::string suffix = fname.substr(prefix_dot.size());
        std::smatch m;
        if (std::regex_match(suffix, m, suffix_pattern)) {
            int k = std::stoi(m[2].str());
            if (filter_k > 0 && k != filter_k) continue;
            VolumeFiles vf;
            vf.volume_index = static_cast<uint16_t>(std::stoi(m[1].str()));
            vf.k = k;
            std::string base = parent_dir + "/" + db_name + "." + m[1].str() + "." + m[2].str() + "mer";
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

// Pre-opened volume data shared across threads (read-only after init).
struct VolumeData {
    KixReader kix;
    KpxReader kpx;
    KsxReader ksx;
    OidFilter filter;
    uint16_t volume_index;
};

// A (query, volume) search job.
struct SearchJob {
    size_t query_idx;
    size_t volume_idx;
};

int main(int argc, char* argv[]) {
    CliParser cli(argc, argv);

    if (cli.has("--version")) {
        std::fprintf(stderr, "ikafssnsearch %s\n", IKAFSSN_VERSION);
        return 0;
    }

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

    std::string ix_prefix = cli.get_string("-ix");
    std::string query_path = cli.get_string("-query");
    int filter_k = cli.get_int("-k", 0);
    std::string output_path = cli.get_string("-o");
    bool verbose = cli.has("-v") || cli.has("--verbose");

    int num_threads = cli.get_int("-threads", 0);
    if (num_threads <= 0) {
        num_threads = static_cast<int>(std::thread::hardware_concurrency());
        if (num_threads <= 0) num_threads = 1;
    }

    Logger logger(verbose ? Logger::kDebug : Logger::kInfo);

    // Search config
    SearchConfig config;
    double max_freq_raw = cli.get_double("-max_freq", 0.5);
    config.stage1.stage1_topn = static_cast<uint32_t>(cli.get_int("-stage1_topn", 0));
    config.stage1.stage1_score_type = static_cast<uint8_t>(cli.get_int("-stage1_score", 1));
    config.stage2.max_gap = static_cast<uint32_t>(cli.get_int("-max_gap", 100));
    config.stage2.min_diag_hits = static_cast<uint32_t>(cli.get_int("-min_diag_hits", 2));
    config.stage2.min_score = static_cast<uint32_t>(cli.get_int("-min_score", 1));
    config.num_results = static_cast<uint32_t>(cli.get_int("-num_results", 0));
    config.mode = static_cast<uint8_t>(cli.get_int("-mode", 2));
    config.sort_score = static_cast<uint8_t>(cli.get_int("-sort_score", 2));
    config.strand = static_cast<int8_t>(cli.get_int("-strand", 2));
    if (config.strand != -1 && config.strand != 1 && config.strand != 2) {
        std::fprintf(stderr, "Error: -strand must be -1, 1, or 2\n");
        return 1;
    }

    // Parse -min_stage1_score as double to support fractional values
    {
        double min_s1 = cli.get_double("-min_stage1_score", 0.5);
        if (min_s1 > 0 && min_s1 < 1.0) {
            config.min_stage1_score_frac = min_s1;
            // Leave stage1.min_stage1_score at default (will be overridden per-query)
        } else {
            config.stage1.min_stage1_score = static_cast<uint32_t>(min_s1);
        }
    }

    // Mode 1: force sort_score=1
    if (config.mode == 1) {
        config.sort_score = 1;

        // Consistency check: fractional + -min_score in mode 1
        if (config.min_stage1_score_frac > 0 && cli.has("-min_score")) {
            std::fprintf(stderr, "Error: -min_score and fractional -min_stage1_score cannot both be specified in -mode 1\n");
            return 1;
        }

        // Consistency check: -min_score and -min_stage1_score in mode 1
        bool has_min_score = cli.has("-min_score");
        bool has_min_stage1_score = cli.has("-min_stage1_score");
        if (config.min_stage1_score_frac == 0) {
            if (has_min_score && has_min_stage1_score &&
                config.stage2.min_score != config.stage1.min_stage1_score) {
                std::fprintf(stderr, "Error: -min_score and -min_stage1_score must be the same in -mode 1\n");
                return 1;
            }
            if (has_min_score && !has_min_stage1_score) {
                config.stage1.min_stage1_score = config.stage2.min_score;
            }
            if (!has_min_score && has_min_stage1_score) {
                config.stage2.min_score = config.stage1.min_stage1_score;
            }
        }
    }

    int accept_qdegen = cli.get_int("-accept_qdegen", 0);

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
    auto vol_files = discover_volumes(ix_prefix, filter_k);
    if (vol_files.empty()) {
        if (filter_k > 0) {
            std::fprintf(stderr, "Error: no index files found for prefix %s with k=%d\n",
                         ix_prefix.c_str(), filter_k);
        } else {
            std::fprintf(stderr, "Error: no index files found for prefix %s\n", ix_prefix.c_str());
        }
        return 1;
    }

    // Determine k: if -k was specified, use it; otherwise require exactly one k value
    int k;
    if (filter_k > 0) {
        k = filter_k;
    } else {
        std::set<int> k_values;
        for (const auto& vf : vol_files) k_values.insert(vf.k);
        if (k_values.size() == 1) {
            k = *k_values.begin();
        } else {
            std::fprintf(stderr, "Error: multiple k-mer sizes found (");
            bool first = true;
            for (int kv : k_values) {
                if (!first) std::fprintf(stderr, ", ");
                std::fprintf(stderr, "%d", kv);
                first = false;
            }
            std::fprintf(stderr, "); specify -k to select one\n");
            return 1;
        }
    }
    logger.info("Found %zu volume(s), k=%d, threads=%d", vol_files.size(), k, num_threads);

    // Read query FASTA
    auto queries = read_fasta(query_path);
    if (queries.empty()) {
        std::fprintf(stderr, "Error: no query sequences found\n");
        return 1;
    }
    logger.info("Read %zu query sequence(s)", queries.size());

    // Check for degenerate bases if accept_qdegen == 0
    std::vector<bool> query_skipped(queries.size(), false);
    bool has_skipped = false;
    if (accept_qdegen == 0) {
        for (size_t qi = 0; qi < queries.size(); qi++) {
            if (contains_degenerate_base(queries[qi].sequence)) {
                query_skipped[qi] = true;
                has_skipped = true;
                std::fprintf(stderr, "Warning: query '%s' contains degenerate bases, skipping\n",
                             queries[qi].id.c_str());
            }
        }
    }

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

    // Pre-open all volumes (mmap kix, optionally kpx, load ksx)
    std::vector<VolumeData> vol_data(vol_files.size());
    for (size_t vi = 0; vi < vol_files.size(); vi++) {
        const auto& vf = vol_files[vi];
        if (!vol_data[vi].kix.open(vf.kix_path)) {
            std::fprintf(stderr, "Error: cannot open %s\n", vf.kix_path.c_str());
            return 1;
        }
        if (config.mode != 1) {
            if (!vol_data[vi].kpx.open(vf.kpx_path)) {
                std::fprintf(stderr, "Error: cannot open %s\n", vf.kpx_path.c_str());
                return 1;
            }
        }
        if (!vol_data[vi].ksx.open(vf.ksx_path)) {
            std::fprintf(stderr, "Error: cannot open %s\n", vf.ksx_path.c_str());
            return 1;
        }
        vol_data[vi].volume_index = vf.volume_index;

        // Build per-volume OID filter
        if (filter_mode != OidFilterMode::kNone) {
            vol_data[vi].filter.build(seqidlist, vol_data[vi].ksx, filter_mode);
        }
    }

    // Open shared .khx (non-fatal if missing)
    KhxReader shared_khx;
    {
        std::filesystem::path prefix_path(ix_prefix);
        std::string parent_dir = prefix_path.parent_path().string();
        std::string db_name = prefix_path.filename().string();
        if (parent_dir.empty()) parent_dir = ".";
        char kk_str[8];
        std::snprintf(kk_str, sizeof(kk_str), "%02d", k);
        std::string khx_path = parent_dir + "/" + db_name + "." + kk_str + "mer.khx";
        shared_khx.open(khx_path); // non-fatal
    }

    // Resolve -max_freq: fraction -> absolute threshold using total NSEQ
    if (max_freq_raw > 0 && max_freq_raw < 1.0) {
        uint64_t total_nseq = 0;
        for (const auto& vd : vol_data) total_nseq += vd.ksx.num_sequences();
        config.stage1.max_freq = static_cast<uint32_t>(
            std::ceil(max_freq_raw * total_nseq));
        if (config.stage1.max_freq == 0) config.stage1.max_freq = 1;
        logger.info("-max_freq=%.6g (fraction) -> threshold=%u (total_nseq=%lu)",
                    max_freq_raw, config.stage1.max_freq,
                    static_cast<unsigned long>(total_nseq));
    } else {
        config.stage1.max_freq = static_cast<uint32_t>(max_freq_raw);
    }

    // Build job list: (query, volume) pairs, skipping degenerate queries
    std::vector<SearchJob> jobs;
    jobs.reserve(queries.size() * vol_data.size());
    for (size_t qi = 0; qi < queries.size(); qi++) {
        if (query_skipped[qi]) continue;
        for (size_t vi = 0; vi < vol_data.size(); vi++) {
            jobs.push_back({qi, vi});
        }
    }

    // Collect results with mutex protection
    std::vector<OutputHit> all_hits;
    std::mutex hits_mutex;

    logger.info("Launching %zu search job(s)...", jobs.size());

    // Execute search jobs in parallel
    tbb::task_arena arena(num_threads);
    arena.execute([&] {
        tbb::parallel_for_each(jobs.begin(), jobs.end(),
            [&](const SearchJob& job) {
                const auto& query = queries[job.query_idx];
                const auto& vd = vol_data[job.volume_idx];

                const KhxReader* khx_ptr = shared_khx.is_open() ? &shared_khx : nullptr;

                SearchResult sr;
                if (k < K_TYPE_THRESHOLD) {
                    sr = search_volume<uint16_t>(
                        query.id, query.sequence, k,
                        vd.kix, vd.kpx, vd.ksx, vd.filter, config, khx_ptr);
                } else {
                    sr = search_volume<uint32_t>(
                        query.id, query.sequence, k,
                        vd.kix, vd.kpx, vd.ksx, vd.filter, config, khx_ptr);
                }

                // Convert to OutputHit and collect
                if (!sr.hits.empty()) {
                    std::vector<OutputHit> local_hits;
                    local_hits.reserve(sr.hits.size());
                    for (const auto& cr : sr.hits) {
                        OutputHit oh;
                        oh.query_id = sr.query_id;
                        oh.accession = std::string(vd.ksx.accession(cr.seq_id));
                        oh.strand = cr.is_reverse ? '-' : '+';
                        oh.q_start = cr.q_start;
                        oh.q_end = cr.q_end;
                        oh.s_start = cr.s_start;
                        oh.s_end = cr.s_end;
                        oh.score = cr.score;
                        oh.stage1_score = cr.stage1_score;
                        oh.volume = vd.volume_index;
                        local_hits.push_back(oh);
                    }

                    std::lock_guard<std::mutex> lock(hits_mutex);
                    all_hits.insert(all_hits.end(), local_hits.begin(), local_hits.end());
                }
            });
    });

    // Sort and truncate final results across volumes (per query)
    if (config.num_results > 0) {
        // Sort by (query_id, sort_score desc)
        if (config.sort_score == 1) {
            std::sort(all_hits.begin(), all_hits.end(),
                      [](const OutputHit& a, const OutputHit& b) {
                          if (a.query_id != b.query_id) return a.query_id < b.query_id;
                          return a.stage1_score > b.stage1_score;
                      });
        } else {
            std::sort(all_hits.begin(), all_hits.end(),
                      [](const OutputHit& a, const OutputHit& b) {
                          if (a.query_id != b.query_id) return a.query_id < b.query_id;
                          return a.score > b.score;
                      });
        }

        // Truncate per query
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
    // num_results == 0: unlimited, skip sort and truncation

    // Write output
    if (output_path.empty()) {
        write_results(std::cout, all_hits, outfmt,
                      config.mode, config.stage1.stage1_score_type);
    } else {
        std::ofstream out(output_path);
        if (!out.is_open()) {
            std::fprintf(stderr, "Error: cannot open output file %s\n", output_path.c_str());
            return 1;
        }
        write_results(out, all_hits, outfmt,
                      config.mode, config.stage1.stage1_score_type);
    }

    logger.info("Done. %zu hit(s) reported.", all_hits.size());
    return has_skipped ? 2 : 0;
}
