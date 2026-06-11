#include "core/config.hpp"
#include "core/spaced_seed.hpp"
#include "core/types.hpp"
#include "core/version.hpp"
#include "util/common_init.hpp"
#include "util/simd_dispatch.hpp"
#include "index/kix_reader.hpp"
#include "index/kpx_reader.hpp"
#include "index/ksx_reader.hpp"
#include "index/khx_reader.hpp"
#include "io/blastdb_reader.hpp"
#include "io/volume_discovery.hpp"
#include "util/cli_parser.hpp"
#include "util/cli_validators.hpp"
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
    print_version_header("ikafssninfo");
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
        "Local mode options (index-variant filter; each unset field is a wildcard):\n"
        "  -db <path>               BLAST DB prefix (default: auto-detect from -ix)\n"
        "  -k <int>                 K-mer length (default: any)\n"
        "  -t <int>                 Template length: 0=contiguous, 13/15/16/18/21\n"
        "                           (default: any; -t 0 selects the contiguous index only)\n"
        "  -template_type <string>  coding, optimal, contiguous, both (default: any; when -t\n"
        "                           is a non-zero value and this is unset, both=coding+optimal)\n"
        "  -min_seq_length <int>    Filter by min_seq_length (default: any)\n"
        "  -min_length_split <int>  Filter by min_length_split (default: any)\n"
        "  -overlap_length <int>    Filter by overlap_length (default: any)\n"
        "  -max_freq_build <int>    Filter by max_freq_build (default: any)\n"
        "  -max_degen_expand <int>  Filter by max_degen_expand (default: any)\n"
        "  -stats <0|1>             Compute k-mer frequency distribution (default: 0). With 1\n"
        "                           this scans the entire posting file and is slow on large\n"
        "                           indexes.\n"
#ifdef IKAFSSN_ENABLE_HTTP
        "\n"
        "Remote HTTP authentication:\n"
        "  -user <user:password>    Credentials (curl-style)\n"
        "  -http_user <USER>        Username (wget-style)\n"
        "  -http_password <PASS>    Password (used with -http_user)\n"
        "  -netrc_file <path>       .netrc file for credentials\n"
#endif
        "\n"
        "Options:\n"
        "  -v, --verbose            Verbose output\n"
        "  -h, --help               Show this help\n",
        prog);
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
    // Per-volume per-kmer posting counts (populated only with -stats)
    std::vector<uint32_t> counts;
    // two-stage KSX surface area.  num_sequences above is
    // the fragment count (= internal SeqId count); num_parents records
    // the parent BLAST OID count behind those fragments.  The two are
    // equal in degenerate (min_length_split == 0) indices and diverge
    // once fragment splitting is enabled.
    uint32_t num_parents = 0;
    // Fragment length distribution (frag_end - frag_start + 1 per
    // internal SeqId).  Empty unless num_sequences > 0.
    uint32_t frag_len_min  = 0;
    uint32_t frag_len_max  = 0;
    double   frag_len_mean = 0.0;
    double   frag_len_median = 0.0;
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

// A one-line descriptor of an index variant, using the same parameters that
// the file name encodes.  Two volumes belong to the same variant iff every
// one of these fields matches.
static std::string variant_label(const DiscoveredVolume& v) {
    std::string s = "k=" + std::to_string(v.k);
    s += " t=" + std::to_string(static_cast<unsigned>(v.t));
    s += " " + template_type_to_string(static_cast<TemplateType>(v.template_type));
    s += " minlen=" + std::to_string(v.min_seq_length);
    s += " minsplit=" + std::to_string(v.min_length_split);
    s += " ovllen=" + std::to_string(v.overlap_length);
    s += " maxfreq=" + std::to_string(static_cast<unsigned long long>(v.max_freq_build));
    s += " maxexpand=" + std::to_string(v.max_degen_expand);
    return s;
}

static bool same_variant(const DiscoveredVolume& a, const DiscoveredVolume& b) {
    return a.k == b.k && a.t == b.t && a.template_type == b.template_type &&
           a.min_seq_length == b.min_seq_length &&
           a.min_length_split == b.min_length_split &&
           a.overlap_length == b.overlap_length &&
           a.max_freq_build == b.max_freq_build &&
           a.max_degen_expand == b.max_degen_expand;
}

static bool variant_less(const DiscoveredVolume& a, const DiscoveredVolume& b) {
    if (a.k != b.k) return a.k < b.k;
    if (a.t != b.t) return a.t < b.t;
    if (a.template_type != b.template_type) return a.template_type < b.template_type;
    if (a.min_seq_length != b.min_seq_length) return a.min_seq_length < b.min_seq_length;
    if (a.min_length_split != b.min_length_split) return a.min_length_split < b.min_length_split;
    if (a.overlap_length != b.overlap_length) return a.overlap_length < b.overlap_length;
    if (a.max_freq_build != b.max_freq_build) return a.max_freq_build < b.max_freq_build;
    if (a.max_degen_expand != b.max_degen_expand) return a.max_degen_expand < b.max_degen_expand;
    return a.volume_index < b.volume_index;
}

// Print the full report for one index variant (all volumes share the same
// encoded parameters).  Statistics are never combined across variants.
static int report_variant(const std::vector<DiscoveredVolume>& vol_files,
                          const std::string& ix_prefix, bool stats_enabled) {
    int k = vol_files[0].k;
    uint8_t vol_t = vol_files[0].t;
    uint8_t vol_template_type = vol_files[0].template_type;
    const uint32_t hdr_min_seq_length   = vol_files[0].min_seq_length;
    const uint32_t hdr_min_length_split = vol_files[0].min_length_split;
    const uint32_t hdr_overlap_length   = vol_files[0].overlap_length;
    uint32_t tbl_size = 0;
    {
        KixReader kix0;
        if (!kix0.open(vol_files[0].kix_path)) {
            std::fprintf(stderr, "Error: cannot open %s\n", vol_files[0].kix_path.c_str());
            return 1;
        }
        tbl_size = kix0.table_size();
        kix0.close();
    }

    // Read all volumes
    std::vector<VolumeStats> vol_stats;
    vol_stats.reserve(vol_files.size());

    uint64_t total_sequences = 0;
    uint64_t total_postings = 0;
    uint64_t total_kix_size = 0;
    uint64_t total_kpx_size = 0;
    uint64_t total_ksx_size = 0;

    // Aggregated counts across all volumes (for frequency distribution).
    // Populated only when -stats is enabled; left empty otherwise.
    std::vector<uint64_t> aggregated_counts;
    if (stats_enabled) aggregated_counts.assign(tbl_size, 0);

    // aggregated parent OID count + fragment-length
    // distribution across all volumes.  Per-volume slices are emitted
    // first (under each "Volume N:" block) and the cross-volume rollup
    // is printed in the overall statistics section.
    uint64_t total_parents = 0;
    std::vector<uint32_t> all_frag_lengths;

    for (const auto& vf : vol_files) {
        KixReader kix;
        if (!kix.open(vf.kix_path)) {
            std::fprintf(stderr, "Error: cannot open %s\n", vf.kix_path.c_str());
            return 1;
        }

        VolumeStats vs;
        vs.volume_index = vf.volume_index;
        vs.num_sequences = kix.num_sequences();
        vs.total_postings = kix.total_distinct_postings();
        vs.kix_size = file_size(vf.kix_path);
        vs.kpx_size = vf.has_kpx ? file_size(vf.kpx_path) : 0;
        vs.ksx_size = file_size(vf.ksx_path);
        vs.has_kpx = vf.has_kpx;

        // Read per-kmer counts for frequency analysis. This scans the entire
        // posting file, so it runs only when -stats is enabled.
        if (stats_enabled) {
            vs.counts = kix.bulk_count_postings();
            for (uint32_t i = 0; i < tbl_size; i++) {
                aggregated_counts[i] += vs.counts[i];
            }
        }

        // open the volume's .ksx to surface the parent /
        // fragment split and the per-volume fragment-length distribution.
        // The KSX reader owns the parent table and per-fragment
        // (start, end) arrays even when the index is degenerate (1 parent
        // == 1 fragment).
        {
            KsxReader ksx;
            if (!ksx.open(vf.ksx_path)) {
                std::fprintf(stderr, "Error: cannot open %s\n", vf.ksx_path.c_str());
                return 1;
            }
            vs.num_parents = ksx.num_parents();

            const uint32_t nseq = ksx.num_sequences();
            if (nseq > 0) {
                std::vector<uint32_t> frag_lens(nseq);
                uint64_t sum_lens = 0;
                vs.frag_len_min = UINT32_MAX;
                vs.frag_len_max = 0;
                for (uint32_t sid = 0; sid < nseq; sid++) {
                    const uint32_t flen =
                        ksx.fragment_end(sid) - ksx.fragment_start(sid) + 1u;
                    frag_lens[sid] = flen;
                    sum_lens += flen;
                    if (flen < vs.frag_len_min) vs.frag_len_min = flen;
                    if (flen > vs.frag_len_max) vs.frag_len_max = flen;
                }
                vs.frag_len_mean = static_cast<double>(sum_lens) /
                                   static_cast<double>(nseq);

                std::vector<uint32_t> sorted_lens = frag_lens;
                std::sort(sorted_lens.begin(), sorted_lens.end());
                const size_t mid = sorted_lens.size() / 2;
                if (sorted_lens.size() % 2 == 0 && !sorted_lens.empty()) {
                    vs.frag_len_median =
                        (static_cast<double>(sorted_lens[mid - 1]) +
                         static_cast<double>(sorted_lens[mid])) / 2.0;
                } else {
                    vs.frag_len_median =
                        static_cast<double>(sorted_lens[mid]);
                }

                all_frag_lengths.insert(all_frag_lengths.end(),
                                        frag_lens.begin(), frag_lens.end());
            }
            ksx.close();
        }
        total_parents += vs.num_parents;

        total_sequences += vs.num_sequences;
        total_postings += vs.total_postings;
        total_kix_size += vs.kix_size;
        total_kpx_size += vs.kpx_size;
        total_ksx_size += vs.ksx_size;

        vol_stats.push_back(std::move(vs));
        kix.close();
    }

    // --- Print index variant information ---
    std::printf("=== Index variant: %s ===\n\n", variant_label(vol_files[0]).c_str());
    std::printf("K-mer length (k):  %d\n", k);
    std::printf("K-mer integer type: %s\n", kmer_type_for(k, vol_t) == 0 ? "uint16" : "uint32");
    if (vol_t > 0) {
        std::printf("Template length:   %u\n", static_cast<unsigned>(vol_t));
        std::printf("Template type:     %s\n",
                     template_type_to_string(static_cast<TemplateType>(vol_template_type)).c_str());
    }
    std::printf("Table size (4^k):  %lu\n", static_cast<unsigned long>(tbl_size));
    std::printf("min_seq_length:    %u\n", hdr_min_seq_length);
    std::printf("min_length_split:  %u\n", hdr_min_length_split);
    std::printf("overlap_length:    %u\n", hdr_overlap_length);
    std::printf("max_freq_build:    %llu\n",
                static_cast<unsigned long long>(vol_files[0].max_freq_build));
    std::printf("max_degen_expand:  %u\n", vol_files[0].max_degen_expand);
    std::printf("Number of volumes: %zu\n\n", vol_stats.size());

    // Per-volume info
    std::printf("--- Per-Volume Statistics ---\n\n");
    for (const auto& vs : vol_stats) {
        std::printf("Volume %u:\n", vs.volume_index);
        // Sequences == fragment count (= internal SeqId
        // count); Parents == parent BLAST OID count.  In a degenerate
        // (min_length_split == 0) index the two are equal.
        std::printf("  Parents:         %u\n", vs.num_parents);
        std::printf("  Fragments:       %u\n", vs.num_sequences);
        if (vs.num_sequences > 0 && vs.num_parents > 0) {
            std::printf("  Frag/Parent:     %.2f\n",
                        static_cast<double>(vs.num_sequences) /
                        static_cast<double>(vs.num_parents));
        }
        if (vs.num_sequences > 0) {
            std::printf("  Fragment length: min=%u median=%.1f mean=%.1f max=%u\n",
                        vs.frag_len_min, vs.frag_len_median,
                        vs.frag_len_mean, vs.frag_len_max);
        }
        std::printf("  Total postings:  %lu\n",
                    static_cast<unsigned long>(vs.total_postings));
        std::printf("  File sizes:\n");
        std::printf("    .kix:          %s (%lu bytes)\n",
                    format_size_human(vs.kix_size).c_str(),
                    static_cast<unsigned long>(vs.kix_size));
        if (vs.has_kpx) {
            std::printf("    .kpx:          %s (%lu bytes)\n",
                        format_size_human(vs.kpx_size).c_str(),
                        static_cast<unsigned long>(vs.kpx_size));
        } else {
            std::printf("    .kpx:          (not built)\n");
        }
        std::printf("    .ksx:          %s (%lu bytes)\n",
                    format_size_human(vs.ksx_size).c_str(),
                    static_cast<unsigned long>(vs.ksx_size));
        uint64_t vol_total = vs.kix_size + vs.kpx_size + vs.ksx_size;
        std::printf("    Total:         %s (%lu bytes)\n",
                    format_size_human(vol_total).c_str(),
                    static_cast<unsigned long>(vol_total));

        if (stats_enabled) {
            FrequencyStats fs = compute_frequency_stats(vs.counts);
            print_frequency_stats(fs);
        }
        std::printf("\n");
    }

    // Try to open the shared .khx.
    auto prefix_parts = parse_index_prefix(ix_prefix);
    const auto& v0 = vol_files[0];
    std::string khx_path = khx_path_for(prefix_parts.parent_dir, prefix_parts.db,
                                         k, vol_t, vol_template_type,
                                         v0.min_seq_length, v0.min_length_split,
                                         v0.overlap_length, v0.max_freq_build,
                                         v0.max_degen_expand);
    KhxReader shared_khx;
    bool has_khx = shared_khx.open(khx_path);
    uint64_t khx_size = has_khx ? file_size(khx_path) : 0;
    uint64_t khx_excluded = has_khx ? shared_khx.count_excluded() : 0;

    if (has_khx) {
        std::printf("--- Shared .khx ---\n\n");
        std::printf("  Path:            %s\n", khx_path.c_str());
        std::printf("  Size:            %s (%lu bytes)\n",
                    format_size_human(khx_size).c_str(),
                    static_cast<unsigned long>(khx_size));
        std::printf("  Excluded k-mers: %lu\n\n",
                    static_cast<unsigned long>(khx_excluded));
    }

    // Overall statistics
    std::printf("--- Overall Statistics ---\n\n");
    // cross-volume parent / fragment rollup.  total_parents
    // and total_sequences (= total fragments) are summed here from the per-
    // volume KSX scans above; the aggregated fragment-length distribution
    // stitches together all per-volume frag_lens vectors.
    std::printf("Total parents:     %lu\n", static_cast<unsigned long>(total_parents));
    std::printf("Total fragments:   %lu\n", static_cast<unsigned long>(total_sequences));
    if (total_parents > 0 && total_sequences > 0) {
        std::printf("Frag/Parent:       %.2f\n",
                    static_cast<double>(total_sequences) /
                    static_cast<double>(total_parents));
    }
    if (!all_frag_lengths.empty()) {
        uint64_t sum_lens = 0;
        uint32_t agg_min = UINT32_MAX;
        uint32_t agg_max = 0;
        for (uint32_t v : all_frag_lengths) {
            sum_lens += v;
            if (v < agg_min) agg_min = v;
            if (v > agg_max) agg_max = v;
        }
        const double agg_mean = static_cast<double>(sum_lens) /
                                static_cast<double>(all_frag_lengths.size());
        std::vector<uint32_t> sorted = all_frag_lengths;
        std::sort(sorted.begin(), sorted.end());
        const size_t mid = sorted.size() / 2;
        const double agg_median =
            (sorted.size() % 2 == 0)
                ? (static_cast<double>(sorted[mid - 1]) +
                   static_cast<double>(sorted[mid])) / 2.0
                : static_cast<double>(sorted[mid]);
        std::printf("Fragment length:   min=%u median=%.1f mean=%.1f max=%u\n",
                    agg_min, agg_median, agg_mean, agg_max);
    }
    std::printf("Total postings:    %lu\n", static_cast<unsigned long>(total_postings));
    uint64_t total_index_size = total_kix_size + total_kpx_size + total_ksx_size + khx_size;
    std::printf("Total index size:  %s (%lu bytes)\n",
                format_size_human(total_index_size).c_str(),
                static_cast<unsigned long>(total_index_size));
    std::printf("  .kix total:      %s\n", format_size_human(total_kix_size).c_str());
    if (total_kpx_size > 0) {
        std::printf("  .kpx total:      %s\n", format_size_human(total_kpx_size).c_str());
    } else {
        std::printf("  .kpx total:      (not built)\n");
    }
    std::printf("  .ksx total:      %s\n", format_size_human(total_ksx_size).c_str());
    if (has_khx) {
        std::printf("  .khx:            %s\n", format_size_human(khx_size).c_str());
    }

    // Compression ratio: compare delta-compressed posting file size vs uncompressed
    // Uncompressed ID posting: total_postings * sizeof(uint32_t) = 4 bytes each
    // Uncompressed pos posting: total_postings * sizeof(uint32_t) = 4 bytes each
    // Total uncompressed posting data: total_postings * (4 or 8) depending on .kpx presence
    bool has_any_kpx = (total_kpx_size > 0);
    uint64_t bytes_per_posting = has_any_kpx ? 8 : 4;  // id + pos, or id only
    uint64_t uncompressed_posting_size = total_postings * bytes_per_posting;
    if (uncompressed_posting_size > 0) {
        // Approximate compressed posting file size: file sizes minus table overhead per volume
        // Per volume: kix header (64) + offsets (8*(4^k+1)) + posting file
        //             kpx (if present): header (32) + pos_offsets (8*4^k) + posting file
        uint64_t table_overhead_per_vol = 64 + (static_cast<uint64_t>(tbl_size) + 1) * 8; // kix
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
        std::printf("  Uncompressed posting file size: %s\n",
                    format_size_human(uncompressed_posting_size).c_str());
        std::printf("  Compressed posting file size:   %s\n",
                    format_size_human(compressed_posting_size).c_str());
        std::printf("  Compression ratio:         %.3f (%.1f%% of original)\n",
                    ratio, ratio * 100.0);
    }

    // Aggregated k-mer frequency distribution
    if (stats_enabled) {
        // Convert aggregated_counts (uint64_t) to uint32_t for stats computation
        // (capped at UINT32_MAX for individual counts, which should not happen in practice)
        std::vector<uint32_t> agg_u32(tbl_size);
        for (uint32_t i = 0; i < tbl_size; i++) {
            agg_u32[i] = (aggregated_counts[i] > UINT32_MAX)
                ? UINT32_MAX
                : static_cast<uint32_t>(aggregated_counts[i]);
        }
        std::printf("\n--- Aggregated K-mer Frequency Distribution ---\n\n");
        FrequencyStats fs = compute_frequency_stats(agg_u32);
        print_frequency_stats(fs);
    }

    std::printf("\n");
    return 0;
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
        if (cli.has("-user") && cli.has("-http_user")) {
            std::fprintf(stderr, "Error: -user and -http_user are mutually exclusive\n");
            return 1;
        }
        if (cli.has("-user")) {
            auth.userpwd = cli.get_string("-user");
        } else if (cli.has("-http_user")) {
            std::string user = cli.get_string("-http_user");
            std::string pass = cli.get_string("-http_password", "");
            auth.userpwd = user + ":" + pass;
        }
        if (cli.has("-netrc_file")) {
            auth.netrc_file = cli.get_string("-netrc_file");
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

    if (cli.has("-h") || cli.has("-help")) {
        print_usage(argv[0]);
        return 0;
    }

    Logger logger = make_logger(cli);
    init_simd_dispatch(&logger);

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
        // -stats scans local posting files and cannot apply to a remote server.
        if (cli.get_int("-stats", 0) != 0) {
            std::fprintf(stderr, "Error: -stats 1 is not supported in remote mode\n");
            return 1;
        }
        // Index-variant filters select among local index files; a remote server
        // exposes a fixed index, so these cannot be applied remotely.
        static const char* const kLocalFilterOpts[] = {
            "-k", "-t", "-template_type", "-min_seq_length", "-min_length_split",
            "-overlap_length", "-max_freq_build", "-max_degen_expand",
        };
        for (const char* opt : kLocalFilterOpts) {
            if (cli.has(opt)) {
                std::fprintf(stderr,
                    "Error: %s is not supported in remote mode\n", opt);
                return 1;
            }
        }
        return run_remote_info(cli, verbose);
    }

    // Local mode
    std::string ix_prefix = cli.get_string("-ix");
    std::string db_path = cli.get_string("-db");

    // Index-variant filter.  Each field is a wildcard unless its option is
    // given; a given field restricts the reported variants to an exact match.
    // -t 0 selects the contiguous index only; a non-zero -t with no
    // -template_type matches both spaced types (coding + optimal).  con/cod/opt
    // (and any other parameter combination) are reported as separate variants,
    // never aggregated into one set of statistics.
    const bool k_set         = cli.has("-k");
    const bool t_set         = cli.has("-t");
    const bool tt_set        = cli.has("-template_type");
    const bool minlen_set    = cli.has("-min_seq_length");
    const bool minsplit_set  = cli.has("-min_length_split");
    const bool ovllen_set    = cli.has("-overlap_length");
    const bool maxfreq_set   = cli.has("-max_freq_build");
    const bool maxexpand_set = cli.has("-max_degen_expand");

    const int filter_k = cli.get_int("-k", 0);

    uint8_t filter_t = 0;
    if (t_set) {
        std::string err;
        if (!parse_spaced_seed_t(cli, filter_t, err)) {
            std::fprintf(stderr, "%s\n", err.c_str());
            return 1;
        }
    }

    TemplateType filter_tt = TemplateType::kContiguous;
    if (tt_set) {
        const std::string s = cli.get_string("-template_type");
        if (s != "coding" && s != "optimal" && s != "contiguous" &&
            s != "both" && s != "coding_and_optimal") {
            std::fprintf(stderr, "Error: invalid -template_type '%s' "
                         "(coding, optimal, contiguous, both)\n", s.c_str());
            return 1;
        }
        filter_tt = template_type_from_string(s);
    }

    const uint32_t filter_minlen    = static_cast<uint32_t>(cli.get_int("-min_seq_length", 0));
    const uint32_t filter_minsplit  = static_cast<uint32_t>(cli.get_int("-min_length_split", 0));
    const uint32_t filter_ovllen    = static_cast<uint32_t>(cli.get_int("-overlap_length", 0));
    const uint64_t filter_maxfreq   = maxfreq_set
        ? std::strtoull(cli.get_string("-max_freq_build").c_str(), nullptr, 10)
        : 0;
    const uint32_t filter_maxexpand = static_cast<uint32_t>(cli.get_int("-max_degen_expand", 0));

    // -stats 1 enables the k-mer frequency distribution, which scans every
    // volume's entire posting file. The default (0) skips that full scan.
    const bool stats_enabled = cli.get_int("-stats", 0) != 0;

    // -t 0 with an explicit -template_type is contradictory.
    {
        std::string err;
        if (!validate_t_template_type_combo(cli, err)) {
            std::fprintf(stderr, "%s\n", err.c_str());
            return 1;
        }
    }

    // Build the shared variant filter (ikafssninfo semantics: an unset -t is a
    // wildcard, and -max_degen_expand is a filter).
    VariantFilter vfilter;
    vfilter.has_k = k_set;                   vfilter.k = filter_k;
    vfilter.t_is_wildcard = !t_set;          vfilter.t = filter_t;
    if (tt_set) {
        vfilter.has_template_type = true;
        vfilter.template_type = static_cast<uint8_t>(filter_tt);
    }
    vfilter.has_min_seq_length = minlen_set;     vfilter.min_seq_length = filter_minlen;
    vfilter.has_min_length_split = minsplit_set; vfilter.min_length_split = filter_minsplit;
    vfilter.has_overlap_length = ovllen_set;     vfilter.overlap_length = filter_ovllen;
    vfilter.has_max_freq_build = maxfreq_set;    vfilter.max_freq_build = filter_maxfreq;
    vfilter.has_max_degen_expand = maxexpand_set; vfilter.max_degen_expand = filter_maxexpand;

    // Discover every volume under the prefix, then apply the variant filter.
    auto all_vols = discover_volumes(ix_prefix);
    std::vector<DiscoveredVolume> matched;
    for (const auto& v : all_vols) {
        if (variant_matches(vfilter, variant_fields_of(v))) matched.push_back(v);
    }

    if (matched.empty()) {
        std::fprintf(stderr, "Error: no index files found for prefix %s "
                     "matching the given filter\n", ix_prefix.c_str());
        return 1;
    }

    // Group the matched volumes into distinct index variants and report each
    // separately; statistics are never combined across variants.
    std::sort(matched.begin(), matched.end(), variant_less);
    std::vector<std::vector<DiscoveredVolume>> variants;
    for (const auto& v : matched) {
        if (variants.empty() || !same_variant(variants.back().front(), v))
            variants.emplace_back();
        variants.back().push_back(v);
    }

    std::printf("=== ikafssn Index Information ===\n\n");
    std::printf("Index prefix:      %s\n", ix_prefix.c_str());
    std::printf("Matched index variants: %zu\n\n", variants.size());

    for (const auto& group : variants) {
        if (report_variant(group, ix_prefix, stats_enabled) != 0)
            return 1;
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
        std::printf("--- BLAST DB Information ---\n\n");
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
