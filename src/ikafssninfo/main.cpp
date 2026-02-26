#include "core/config.hpp"
#include "core/types.hpp"
#include "core/version.hpp"
#include "util/common_init.hpp"
#include "index/kix_reader.hpp"
#include "index/kpx_reader.hpp"
#include "index/ksx_reader.hpp"
#include "index/khx_reader.hpp"
#include "io/blastdb_reader.hpp"
#include "io/volume_discovery.hpp"
#include "util/cli_parser.hpp"
#include "util/logger.hpp"
#include "util/socket_utils.hpp"
#include "ikafssnclient/socket_client.hpp"
#ifdef IKAFSSN_ENABLE_HTTP
#include "ikafssnclient/http_client.hpp"
#endif
#include "protocol/info_format.hpp"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <numeric>
#include <string>
#include <vector>

using namespace ikafssn;

static void print_usage(const char* prog) {
    std::fprintf(stderr,
        "Usage: %s [options]\n"
        "\n"
        "Required (one of):\n"
        "  -ix <prefix>             Index prefix [local mode]\n"
        "  -socket <path>           UNIX socket to ikafssnserver [remote mode]\n"
        "  -tcp <host>:<port>       TCP address of ikafssnserver [remote mode]\n"
#ifdef IKAFSSN_ENABLE_HTTP
        "  -http <url>              ikafssnhttpd URL [remote mode]\n"
#endif
        "\n"
        "Local mode options:\n"
        "  -db <path>               BLAST DB prefix (default: auto-detect from -ix)\n"
#ifdef IKAFSSN_ENABLE_HTTP
        "\n"
        "Remote HTTP authentication:\n"
        "  --user <user:password>   Credentials (curl-style)\n"
        "  --http-user <USER>       Username (wget-style)\n"
        "  --http-password <PASS>   Password (used with --http-user)\n"
        "  --netrc-file <path>      .netrc file for credentials\n"
#endif
        "\n"
        "Options:\n"
        "  -v, --verbose            Verbose output\n"
        "  -h, --help               Show this help\n",
        prog);
}

static std::string format_size(uint64_t bytes) {
    if (bytes >= uint64_t(1) << 30) {
        return std::to_string(bytes / (uint64_t(1) << 30)) + "."
             + std::to_string((bytes % (uint64_t(1) << 30)) * 10 / (uint64_t(1) << 30))
             + " GiB";
    } else if (bytes >= uint64_t(1) << 20) {
        return std::to_string(bytes / (uint64_t(1) << 20)) + "."
             + std::to_string((bytes % (uint64_t(1) << 20)) * 10 / (uint64_t(1) << 20))
             + " MiB";
    } else if (bytes >= uint64_t(1) << 10) {
        return std::to_string(bytes / (uint64_t(1) << 10)) + "."
             + std::to_string((bytes % (uint64_t(1) << 10)) * 10 / (uint64_t(1) << 10))
             + " KiB";
    }
    return std::to_string(bytes) + " B";
}

static uint64_t file_size(const std::string& path) {
    std::error_code ec;
    auto sz = std::filesystem::file_size(path, ec);
    if (ec) return 0;
    return sz;
}

struct VolumeStats {
    uint16_t volume_index;
    uint32_t num_sequences;
    uint64_t total_postings;
    uint64_t kix_size;
    uint64_t kpx_size;
    uint64_t ksx_size;
    bool has_kpx;
    // Per-volume counts array (for verbose stats)
    std::vector<uint32_t> counts;
};

struct FrequencyStats {
    uint64_t total_entries;     // sum of all counts
    uint64_t non_empty_kmers;  // number of k-mers with count > 0
    uint64_t total_kmers;      // 4^k
    uint32_t min_count;
    uint32_t max_count;
    double mean;
    double median;
    double p25;
    double p75;
    double p95;
    double p99;
};

static FrequencyStats compute_frequency_stats(const std::vector<uint32_t>& counts) {
    FrequencyStats fs = {};
    fs.total_kmers = counts.size();
    fs.total_entries = 0;
    fs.min_count = UINT32_MAX;
    fs.max_count = 0;
    fs.non_empty_kmers = 0;

    for (uint32_t c : counts) {
        fs.total_entries += c;
        if (c > 0) {
            fs.non_empty_kmers++;
            if (c < fs.min_count) fs.min_count = c;
            if (c > fs.max_count) fs.max_count = c;
        }
    }

    if (fs.non_empty_kmers == 0) {
        fs.min_count = 0;
        return fs;
    }

    fs.mean = static_cast<double>(fs.total_entries) / static_cast<double>(fs.total_kmers);

    // Compute percentiles: sort counts
    std::vector<uint32_t> sorted_counts(counts.begin(), counts.end());
    std::sort(sorted_counts.begin(), sorted_counts.end());

    auto percentile = [&](double p) -> double {
        double idx = p * static_cast<double>(sorted_counts.size() - 1);
        size_t lo = static_cast<size_t>(idx);
        size_t hi = lo + 1;
        if (hi >= sorted_counts.size()) hi = sorted_counts.size() - 1;
        double frac = idx - static_cast<double>(lo);
        return static_cast<double>(sorted_counts[lo]) * (1.0 - frac)
             + static_cast<double>(sorted_counts[hi]) * frac;
    };

    fs.median = percentile(0.5);
    fs.p25 = percentile(0.25);
    fs.p75 = percentile(0.75);
    fs.p95 = percentile(0.95);
    fs.p99 = percentile(0.99);

    return fs;
}

static void print_frequency_stats(const FrequencyStats& fs) {
    std::printf("  K-mer frequency distribution:\n");
    std::printf("    Total k-mer slots:     %lu (4^k)\n",
                static_cast<unsigned long>(fs.total_kmers));
    std::printf("    Non-empty k-mers:      %lu (%.1f%%)\n",
                static_cast<unsigned long>(fs.non_empty_kmers),
                fs.total_kmers > 0
                    ? 100.0 * static_cast<double>(fs.non_empty_kmers) / static_cast<double>(fs.total_kmers)
                    : 0.0);
    std::printf("    Min count:             %u\n", fs.min_count);
    std::printf("    Max count:             %u\n", fs.max_count);
    std::printf("    Mean count:            %.2f\n", fs.mean);
    std::printf("    Percentiles:\n");
    std::printf("      25th:                %.1f\n", fs.p25);
    std::printf("      50th (median):       %.1f\n", fs.median);
    std::printf("      75th:                %.1f\n", fs.p75);
    std::printf("      95th:                %.1f\n", fs.p95);
    std::printf("      99th:                %.1f\n", fs.p99);
}

static int run_remote_info(const CliParser& cli, bool verbose) {
    InfoResponse info;

    if (cli.has("-socket") || cli.has("-tcp")) {
        int fd = -1;
        if (cli.has("-socket")) {
            std::string sock_path = cli.get_string("-socket");
            fd = unix_connect(sock_path);
            if (fd < 0) {
                std::fprintf(stderr, "Error: cannot connect to UNIX socket %s\n",
                             sock_path.c_str());
                return 1;
            }
        } else {
            std::string tcp_addr = cli.get_string("-tcp");
            fd = tcp_connect(tcp_addr);
            if (fd < 0) {
                std::fprintf(stderr, "Error: cannot connect to TCP %s\n",
                             tcp_addr.c_str());
                return 1;
            }
        }

        if (!socket_info(fd, info)) {
            std::fprintf(stderr, "Error: info request failed\n");
            close_fd(fd);
            return 1;
        }
        close_fd(fd);
    }
#ifdef IKAFSSN_ENABLE_HTTP
    else if (cli.has("-http")) {
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

        std::string error_msg;
        if (!http_info(cli.get_string("-http"), info, error_msg, auth)) {
            std::fprintf(stderr, "Error: %s\n", error_msg.c_str());
            return 1;
        }
    }
#endif
    else {
        std::fprintf(stderr, "Error: no remote connection specified\n");
        return 1;
    }

    std::printf("%s", format_server_info(info, verbose).c_str());
    return 0;
}

int main(int argc, char* argv[]) {
    CliParser cli(argc, argv);

    if (check_version(cli, "ikafssninfo")) return 0;

    if (cli.has("-h") || cli.has("--help")) {
        print_usage(argv[0]);
        return 0;
    }

    bool has_ix = cli.has("-ix");
    bool has_socket = cli.has("-socket");
    bool has_tcp = cli.has("-tcp");
    bool has_http = false;
#ifdef IKAFSSN_ENABLE_HTTP
    has_http = cli.has("-http");
#endif
    bool verbose = cli.has("-v") || cli.has("--verbose");

    bool has_remote = has_socket || has_tcp || has_http;

    if (!has_ix && !has_remote) {
        print_usage(argv[0]);
        return 1;
    }

    // -ix and remote options are mutually exclusive
    if (has_ix && has_remote) {
        std::fprintf(stderr, "Error: -ix cannot be used with remote options "
                     "(-socket, -tcp"
#ifdef IKAFSSN_ENABLE_HTTP
                     ", -http"
#endif
                     ")\n");
        return 1;
    }

    // Only one remote option at a time
    if ((has_socket ? 1 : 0) + (has_tcp ? 1 : 0) + (has_http ? 1 : 0) > 1) {
        std::fprintf(stderr, "Error: only one remote option "
                     "(-socket, -tcp"
#ifdef IKAFSSN_ENABLE_HTTP
                     ", -http"
#endif
                     ") may be specified\n");
        return 1;
    }

    // Remote mode
    if (has_remote) {
        return run_remote_info(cli, verbose);
    }

    // Local mode
    std::string ix_prefix = cli.get_string("-ix");
    std::string db_path = cli.get_string("-db");

    // Discover volumes
    auto vol_files = discover_volumes(ix_prefix);
    if (vol_files.empty()) {
        std::fprintf(stderr, "Error: no index files found for prefix %s\n", ix_prefix.c_str());
        return 1;
    }

    int k = vol_files[0].k;
    uint64_t tbl_size = table_size(k);

    // Read all volumes
    std::vector<VolumeStats> vol_stats;
    vol_stats.reserve(vol_files.size());

    uint64_t total_sequences = 0;
    uint64_t total_postings = 0;
    uint64_t total_kix_size = 0;
    uint64_t total_kpx_size = 0;
    uint64_t total_ksx_size = 0;

    // Aggregated counts across all volumes (for frequency distribution)
    std::vector<uint64_t> aggregated_counts(tbl_size, 0);

    for (const auto& vf : vol_files) {
        KixReader kix;
        if (!kix.open(vf.kix_path)) {
            std::fprintf(stderr, "Error: cannot open %s\n", vf.kix_path.c_str());
            return 1;
        }

        VolumeStats vs;
        vs.volume_index = vf.volume_index;
        vs.num_sequences = kix.num_sequences();
        vs.total_postings = kix.total_postings();
        vs.kix_size = file_size(vf.kix_path);
        vs.kpx_size = vf.has_kpx ? file_size(vf.kpx_path) : 0;
        vs.ksx_size = file_size(vf.ksx_path);
        vs.has_kpx = vf.has_kpx;

        // Read per-kmer counts for frequency analysis
        const uint32_t* cts = kix.counts();
        vs.counts.assign(cts, cts + tbl_size);

        for (uint64_t i = 0; i < tbl_size; i++) {
            aggregated_counts[i] += cts[i];
        }

        total_sequences += vs.num_sequences;
        total_postings += vs.total_postings;
        total_kix_size += vs.kix_size;
        total_kpx_size += vs.kpx_size;
        total_ksx_size += vs.ksx_size;

        vol_stats.push_back(std::move(vs));
        kix.close();
    }

    // --- Print index information ---
    std::printf("=== ikafssn Index Information ===\n\n");
    std::printf("Index prefix:      %s\n", ix_prefix.c_str());
    std::printf("K-mer length (k):  %d\n", k);
    std::printf("K-mer integer type: %s\n", k < K_TYPE_THRESHOLD ? "uint16" : "uint32");
    std::printf("Table size (4^k):  %lu\n", static_cast<unsigned long>(tbl_size));
    std::printf("Number of volumes: %zu\n\n", vol_stats.size());

    // Per-volume info
    std::printf("--- Per-Volume Statistics ---\n\n");
    for (const auto& vs : vol_stats) {
        std::printf("Volume %u:\n", vs.volume_index);
        std::printf("  Sequences:       %u\n", vs.num_sequences);
        std::printf("  Total postings:  %lu\n",
                    static_cast<unsigned long>(vs.total_postings));
        std::printf("  File sizes:\n");
        std::printf("    .kix:          %s (%lu bytes)\n",
                    format_size(vs.kix_size).c_str(),
                    static_cast<unsigned long>(vs.kix_size));
        if (vs.has_kpx) {
            std::printf("    .kpx:          %s (%lu bytes)\n",
                        format_size(vs.kpx_size).c_str(),
                        static_cast<unsigned long>(vs.kpx_size));
        } else {
            std::printf("    .kpx:          (not built)\n");
        }
        std::printf("    .ksx:          %s (%lu bytes)\n",
                    format_size(vs.ksx_size).c_str(),
                    static_cast<unsigned long>(vs.ksx_size));
        uint64_t vol_total = vs.kix_size + vs.kpx_size + vs.ksx_size;
        std::printf("    Total:         %s (%lu bytes)\n",
                    format_size(vol_total).c_str(),
                    static_cast<unsigned long>(vol_total));

        if (verbose) {
            FrequencyStats fs = compute_frequency_stats(vs.counts);
            print_frequency_stats(fs);
        }
        std::printf("\n");
    }

    // Try to open shared .khx
    auto prefix_parts = parse_index_prefix(ix_prefix);
    std::string khx_path = khx_path_for(prefix_parts.parent_dir, prefix_parts.db_name, k);
    KhxReader shared_khx;
    bool has_khx = shared_khx.open(khx_path);
    uint64_t khx_size = has_khx ? file_size(khx_path) : 0;
    uint64_t khx_excluded = has_khx ? shared_khx.count_excluded() : 0;

    if (has_khx) {
        std::printf("--- Shared .khx ---\n\n");
        std::printf("  Path:            %s\n", khx_path.c_str());
        std::printf("  Size:            %s (%lu bytes)\n",
                    format_size(khx_size).c_str(),
                    static_cast<unsigned long>(khx_size));
        std::printf("  Excluded k-mers: %lu\n\n",
                    static_cast<unsigned long>(khx_excluded));
    }

    // Overall statistics
    std::printf("--- Overall Statistics ---\n\n");
    std::printf("Total sequences:   %lu\n", static_cast<unsigned long>(total_sequences));
    std::printf("Total postings:    %lu\n", static_cast<unsigned long>(total_postings));
    uint64_t total_index_size = total_kix_size + total_kpx_size + total_ksx_size + khx_size;
    std::printf("Total index size:  %s (%lu bytes)\n",
                format_size(total_index_size).c_str(),
                static_cast<unsigned long>(total_index_size));
    std::printf("  .kix total:      %s\n", format_size(total_kix_size).c_str());
    if (total_kpx_size > 0) {
        std::printf("  .kpx total:      %s\n", format_size(total_kpx_size).c_str());
    } else {
        std::printf("  .kpx total:      (not built)\n");
    }
    std::printf("  .ksx total:      %s\n", format_size(total_ksx_size).c_str());
    if (has_khx) {
        std::printf("  .khx:            %s\n", format_size(khx_size).c_str());
    }

    // Compression ratio: compare delta-compressed posting size vs uncompressed
    // Uncompressed ID posting: total_postings * sizeof(uint32_t) = 4 bytes each
    // Uncompressed pos posting: total_postings * sizeof(uint32_t) = 4 bytes each
    // Total uncompressed posting data: total_postings * (4 or 8) depending on .kpx presence
    bool has_any_kpx = (total_kpx_size > 0);
    uint64_t bytes_per_posting = has_any_kpx ? 8 : 4;  // id + pos, or id only
    uint64_t uncompressed_posting_size = total_postings * bytes_per_posting;
    if (uncompressed_posting_size > 0) {
        // Approximate compressed posting size: file sizes minus table overhead per volume
        // Per volume: kix header (64) + offsets (8*4^k) + counts (4*4^k) + posting data
        //             kpx (if present): header (32) + pos_offsets (8*4^k) + posting data
        uint64_t table_overhead_per_vol = 64 + tbl_size * 8 + tbl_size * 4; // kix
        if (has_any_kpx) {
            table_overhead_per_vol += 32 + tbl_size * 8;                     // kpx
        }
        uint64_t total_table_overhead = table_overhead_per_vol * vol_stats.size();
        uint64_t compressed_posting_size = (total_kix_size + total_kpx_size > total_table_overhead)
            ? (total_kix_size + total_kpx_size - total_table_overhead)
            : 0;
        double ratio = static_cast<double>(compressed_posting_size)
                     / static_cast<double>(uncompressed_posting_size);
        std::printf("\nCompression:\n");
        std::printf("  Uncompressed posting size: %s\n",
                    format_size(uncompressed_posting_size).c_str());
        std::printf("  Compressed posting size:   %s\n",
                    format_size(compressed_posting_size).c_str());
        std::printf("  Compression ratio:         %.3f (%.1f%% of original)\n",
                    ratio, ratio * 100.0);
    }

    // Aggregated k-mer frequency distribution
    if (verbose) {
        // Convert aggregated_counts (uint64_t) to uint32_t for stats computation
        // (capped at UINT32_MAX for individual counts, which should not happen in practice)
        std::vector<uint32_t> agg_u32(tbl_size);
        for (uint64_t i = 0; i < tbl_size; i++) {
            agg_u32[i] = (aggregated_counts[i] > UINT32_MAX)
                ? UINT32_MAX
                : static_cast<uint32_t>(aggregated_counts[i]);
        }
        std::printf("\n--- Aggregated K-mer Frequency Distribution ---\n\n");
        FrequencyStats fs = compute_frequency_stats(agg_u32);
        print_frequency_stats(fs);
    }

    // BLAST DB information: use -db if specified, otherwise try -ix path
    std::string effective_db = db_path;
    if (effective_db.empty()) {
        auto vol_paths = BlastDbReader::find_volume_paths(ix_prefix);
        if (!vol_paths.empty()) {
            effective_db = ix_prefix;
        }
    }

    if (!effective_db.empty()) {
        std::printf("\n--- BLAST DB Information ---\n\n");
        std::printf("DB prefix:         %s\n", effective_db.c_str());

        auto vol_paths = BlastDbReader::find_volume_paths(effective_db);
        std::printf("DB volumes:        %zu\n", vol_paths.size());

        // Open the DB to get title and aggregate stats
        BlastDbReader db;
        if (db.open(effective_db)) {
            std::string title = db.get_title();
            uint32_t db_nseqs = db.num_sequences();
            std::printf("DB title:          %s\n", title.c_str());
            std::printf("DB sequences:      %u\n", db_nseqs);

            // Compute total bases
            uint64_t total_bases = 0;
            for (uint32_t oid = 0; oid < db_nseqs; oid++) {
                total_bases += db.seq_length(oid);
            }
            std::printf("DB total bases:    %lu\n",
                        static_cast<unsigned long>(total_bases));

            if (verbose && !vol_paths.empty()) {
                std::printf("\n  DB volume paths:\n");
                for (size_t i = 0; i < vol_paths.size(); i++) {
                    std::printf("    [%zu] %s\n", i, vol_paths[i].c_str());
                }
            }
            db.close();
        } else {
            std::fprintf(stderr, "Warning: could not open BLAST DB '%s'\n",
                         effective_db.c_str());
        }
    }

    std::printf("\n");
    return 0;
}
