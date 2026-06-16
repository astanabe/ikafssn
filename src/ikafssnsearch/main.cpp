#include "core/config.hpp"
#include "core/spaced_seed.hpp"
#include "core/types.hpp"
#include "core/version.hpp"
#include "core/kmer_encoding.hpp"
#include "io/primer_query.hpp"
#include "protocol/messages.hpp"
#include "index/kix_reader.hpp"
#include "index/kpx_reader.hpp"
#include "index/ksx_reader.hpp"
#include "index/khx_reader.hpp"
#include "search/oid_filter.hpp"
#include "search/result_dedup.hpp"
#include "search/search_config.hpp"
#include "search/parallel_search.hpp"
#include "search/search_orchestrator.hpp"
#include "search/preprocess_runner.hpp"
#include "search/query_preprocessor.hpp"
#include "search/stage3_alignment.hpp"
#include "search/width_selection.hpp"
#include "io/fasta_reader.hpp"
#include "io/blastdb_reader.hpp"
#include "io/compressed_stream.hpp"
#include "io/volume_discovery.hpp"
#include "io/seqidlist_reader.hpp"
#include "io/result_writer.hpp"
#include "io/sam_writer.hpp"
#include "util/cli_parser.hpp"
#include "util/cli_validators.hpp"
#include "util/common_init.hpp"
#include "util/simd_dispatch.hpp"
#include "util/context_parser.hpp"
#include "util/logger.hpp"
#include "util/size_parser.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <thread>
#include <vector>

#include <tbb/blocked_range.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>
#include <tbb/task_arena.h>

using namespace ikafssn;

static void print_usage(const char* prog) {
    print_version_header("ikafssnsearch");
    std::fprintf(stderr,
        "Usage: %s [options]\n"
        "\n"
        "Required:\n"
        "  -ix <prefix>             Index prefix (like blastn -db)\n"
        "  -query <path>            Query FASTA file (- for stdin)\n"
        "\n"
        "Primer mode (alternative to -query):\n"
        "  -primer <path>           Primer pair FASTA (even number of sequences; mutually exclusive with -query)\n"
        "  -insert_length <int>     Expected insert length (required with -primer)\n"
        "  -stage1_primer_score <num>  Stage 1 threshold for primer (0<v<=1: fraction, v>=2: absolute; default: 0.5)\n"
        "  -stage2_primer_score_add <int>  Stage 2 threshold addon: max(Lf,Lr) + N (default: 1)\n"
        "\n"
        "Options:\n"
        "  -k <int>                 K-mer size to use (required if multiple k values exist)\n"
        "  -o <path>                Output file (default: stdout)\n"
        "  -nthread <int>           Parallel search threads (default: all cores)\n"
        "  -memory_limit <size>     Search memory budget (default: half of RAM)\n"
        "                           Pins .khx/.ksx metadata; residual caps the\n"
        "                           Stage 3 posting heap.  Accepts K, M, G suffixes\n"
        "  -mode <1|2|3>            1=Stage1, 2=Stage1+2, 3=Stage1+2+3 (default: 1)\n"
        "  -min_query_length <int>  Minimum query length; shorter queries are skipped with a warning.\n"
        "                           Must be >= the index's min_seq_length (default: 64)\n"
        "  -db <path>               BLAST DB path for mode 3 (default: same as -ix)\n"
        "  -stage2_min_score <int>  Minimum chain score (default: 0 = adaptive)\n"
        "                           0 = use resolved Stage 1 threshold\n"
        "  -stage2_max_gap <int>    Chaining diagonal gap tolerance (default: 100)\n"
        "  -stage2_max_lookback <int>  Chaining DP lookback window (default: 64, 0=unlimited)\n"
        "  -stage2_max_nhit_per_subject <int>  Max chains per subject (default: 1, 0=unlimited)\n"
        "  -stage2_max_nhit_per_subject_mode <1|2|3|4>  Per-subject selection mode (default: 3)\n"
        "                           1/2=top-N (no ties), 3/4=top-N + score ties;\n"
        "                           1/3=strands merged per parent, 2/4=strands separate\n"
        "  -stage2_max_nhit_in_total <int>  Stage 2: max chains per query over all volumes,\n"
        "                           by chainscore (default: 0=unlimited)\n"
        "  -stage2_min_nhit_diag <int>  Diagonal filter min hits (default: 1)\n"
        "  -stage1_max_nhit_per_subject <int>  Stage 1: max candidates per parent (default: 1, 0=unlimited)\n"
        "  -stage1_max_nhit_per_subject_mode <1|2|3|4>  Per-subject selection mode (default: 3)\n"
        "                           1/2=top-N (no ties), 3/4=top-N + score ties;\n"
        "                           1/3=strands merged per parent, 2/4=strands separate\n"
        "  -stage1_max_nhit_per_volume <int>  Stage 1: max candidates per (query,volume,strand)\n"
        "                           (default: 0=unlimited)\n"
        "  -stage1_max_nhit_in_total <int>  Stage 1: max candidates per query over all volumes;\n"
        "                           setting this also sets -stage1_max_nhit_per_volume to the same\n"
        "                           value unless that is given explicitly (must be >= it) (default: 0=unlimited)\n"
        "  -stage1_min_score <num>  Stage 1 minimum score; integer or 0<P<1 fraction (default: 0.5)\n"
        "  -seqidlist <path>        Include only listed accessions\n"
        "  -negative_seqidlist <path>  Exclude listed accessions\n"
        "  -strand <-1|1|2>         Strand: 1=plus, -1=minus, 2=both (default: 2)\n"
        "  -accept_qdegen <0|1>     Accept queries with degenerate bases (default: 1)\n"
        "  -context_extend <value>  Context extension for mode 3 (int=bases, decimal=query length multiplier, default: 2.0)\n"
        "  -stage3_traceback <0|1>  Enable traceback in mode 3 (default: 0)\n"
        "  -stage3_gapopen <int>    Gap open penalty for mode 3 (default: 10)\n"
        "  -stage3_gapext <int>     Gap extension penalty for mode 3 (default: 1)\n"
        "  -stage3_min_ppositive <num> Min percent positive filter for mode 3 (default: 0)\n"
        "  -stage3_min_npositive <int> Min positive-scoring positions filter for mode 3 (default: 0)\n"
        "  -stage3_max_nhit_per_subject <int>  Stage 3: max hits per subject, by alnscore (default: 1, 0=unlimited)\n"
        "  -stage3_max_nhit_per_subject_mode <1|2|3|4>  Per-subject selection mode (default: 3)\n"
        "                           1/2=top-N (no ties), 3/4=top-N + score ties;\n"
        "                           1/3=strands merged per parent, 2/4=strands separate\n"
        "  -stage3_max_nhit_in_total <int>  Stage 3: max hits per query over all volumes,\n"
        "                           by alnscore (default: 0=unlimited)\n"
        "  -stage3_score_matrix <name>  Score matrix: degmatch, dnafull, nuc44 (default: degmatch)\n"
        "  -max_degen_expand <int>  Max degenerate expansion per k-mer (default: 16, max: 256, 0/1: disable)\n"
        "  -t <int>                 Template length for spaced seeds (0=contiguous, 13/15/18 for k=8-9, 16/18/21 for k=11-12; default: 0)\n"
        "  -template_type <string>  Template type: coding, optimal, both (default: both; invalid with -t 0)\n"
        "  -min_seq_length <int>    Select the index variant with this min_seq_length (default: any)\n"
        "  -min_length_split <int>  Select the index variant with this min_length_split (default: any)\n"
        "  -overlap_length <int>    Select the index variant with this overlap_length (default: any)\n"
        "  -max_freq_build <int>    Select the index variant with this max_freq_build (default: any)\n"
        "  -output_format <tsv|json|sam|bam>  Output format (default: tsv)\n"
        "  -compression_level <int>    Output compression level (default per codec: gzip=6, bzip2=9, xz=6, zstd=3)\n"
        "  -v, --verbose            Verbose logging\n",
        prog);
}

// Per-volume static metadata: paths, posting/full sizes, ksx reader.  The
// .kix / .kpx readers are opened/closed per-batch inside run_search; only
// .ksx is held open here for accession / seq_length lookups during
// OutputHit conversion (and for OidFilter construction).
struct VolumeData {
    DiscoveredVolume files;
    KsxReader ksx;
    OidFilter filter;
    uint32_t num_sequences = 0;
    uint16_t volume_index = 0;
};

// Per-template-side data. For "both" template_type the search holds two
// contexts (coding + optimal); otherwise a single context is used. Each
// context owns its own volume metadata, KHX bitset, seed masks, and
// preprocessed query data, so the rest of the pipeline iterates uniformly
// over `ctxs` rather than duplicating coding/optimal logic.
template <typename KmerInt>
struct PreprocessedQuery { QueryKmerData<KmerInt> qdata; };

struct TemplateContext {
    TemplateType type = TemplateType::kBoth;
    std::vector<DiscoveredVolume> vol_files;
    std::vector<VolumeData> volumes;
    KhxReader khx;
    const KhxReader* khx_ptr = nullptr;
    std::vector<uint32_t> seed_masks;
    std::vector<PreprocessedQuery<uint16_t>> pp16;
    std::vector<PreprocessedQuery<uint32_t>> pp32;
};

int main(int argc, char* argv[]) {
    CliParser cli(argc, argv);

    if (check_version(cli, "ikafssnsearch")) return 0;

    if (cli.has("-h") || cli.has("-help")) {
        print_usage(argv[0]);
        return 0;
    }

    bool has_query = cli.has("-query");
    bool has_primer = cli.has("-primer");

    if (!cli.has("-ix") || (!has_query && !has_primer)) {
        print_usage(argv[0]);
        return 1;
    }

    if (has_query && has_primer) {
        std::fprintf(stderr, "Error: -query and -primer are mutually exclusive\n");
        return 1;
    }

    // Check mutually exclusive options
    if (cli.has("-seqidlist") && cli.has("-negative_seqidlist")) {
        std::fprintf(stderr, "Error: -seqidlist and -negative_seqidlist are mutually exclusive\n");
        return 1;
    }

    // Primer mode validation
    {
        std::string err;
        if (!validate_primer_mode_options(cli, err)) {
            std::fprintf(stderr, "%s\n", err.c_str());
            return 1;
        }
    }

    // -t 0 with an explicit -template_type is contradictory.
    {
        std::string err;
        if (!validate_t_template_type_combo(cli, err)) {
            std::fprintf(stderr, "%s\n", err.c_str());
            return 1;
        }
    }

    std::string ix_prefix = cli.get_string("-ix");
    std::string query_path = has_query ? cli.get_string("-query") : "";
    std::string primer_path = has_primer ? cli.get_string("-primer") : "";
    int filter_k = cli.get_int("-k", 0);
    std::string output_path = cli.get_string("-o");
    int nthread = resolve_threads(cli);
    Logger logger = make_logger(cli);
    init_simd_dispatch(&logger);

    // Memory limit for madvise budget
    uint64_t memory_limit;
    if (cli.has("-memory_limit")) {
        std::string mem_str = cli.get_string("-memory_limit");
        memory_limit = parse_size_string(mem_str);
        if (memory_limit == 0) {
            std::fprintf(stderr, "Error: invalid -memory_limit '%s'\n", mem_str.c_str());
            return 1;
        }
    } else {
        memory_limit = default_memory_limit();
    }

    // Search config
    SearchConfig config;
    config.stage1.max_nhit_per_subject =
        static_cast<uint32_t>(cli.get_int("-stage1_max_nhit_per_subject", 1));
    if (cli.has("-stage1_max_nhit_per_subject_mode")) {
        int m = cli.get_int("-stage1_max_nhit_per_subject_mode", 3);
        if (m < 1 || m > 4) {
            std::fprintf(stderr,
                "Error: -stage1_max_nhit_per_subject_mode must be 1, 2, 3, or 4\n");
            return 1;
        }
        config.stage1.max_nhit_per_subject_mode = static_cast<uint8_t>(m);
    }
    // M / L per-volume / in-total candidate caps.  Setting only L copies it to M;
    // setting only M leaves L unlimited.  When both are > 0, require L >= M.
    {
        const bool has_m = cli.has("-stage1_max_nhit_per_volume");
        const bool has_l = cli.has("-stage1_max_nhit_in_total");
        uint32_t mv = static_cast<uint32_t>(cli.get_int("-stage1_max_nhit_per_volume", 0));
        uint32_t lv = static_cast<uint32_t>(cli.get_int("-stage1_max_nhit_in_total", 0));
        if (has_l && !has_m) mv = lv;
        if (mv != 0 && lv != 0 && lv < mv) {
            std::fprintf(stderr,
                "Error: -stage1_max_nhit_in_total must be >= -stage1_max_nhit_per_volume\n");
            return 1;
        }
        config.stage1.max_nhit_per_volume = mv;
        config.stage1.max_nhit_in_total = lv;
    }
    config.stage2.max_gap = static_cast<uint32_t>(cli.get_int("-stage2_max_gap", 100));
    config.stage2.chain_max_lookback = static_cast<uint32_t>(cli.get_int("-stage2_max_lookback", 64));
    config.stage2.max_nhit_per_subject = static_cast<uint32_t>(cli.get_int("-stage2_max_nhit_per_subject", 1));
    if (cli.has("-stage2_max_nhit_per_subject_mode")) {
        int m = cli.get_int("-stage2_max_nhit_per_subject_mode", 3);
        if (m < 1 || m > 4) {
            std::fprintf(stderr,
                "Error: -stage2_max_nhit_per_subject_mode must be 1, 2, 3, or 4\n");
            return 1;
        }
        config.stage2.max_nhit_per_subject_mode = static_cast<uint8_t>(m);
    }
    config.stage2.max_nhit_in_total =
        static_cast<uint32_t>(cli.get_int("-stage2_max_nhit_in_total", 0));
    config.stage2.min_nhit_diag = static_cast<uint32_t>(cli.get_int("-stage2_min_nhit_diag", 1));
    config.stage2.min_score = static_cast<uint32_t>(cli.get_int("-stage2_min_score", 0));
    config.mode = static_cast<uint8_t>(cli.get_int("-mode", 1));
    config.strand = static_cast<int8_t>(cli.get_int("-strand", 2));
    if (config.strand != -1 && config.strand != 1 && config.strand != 2) {
        std::fprintf(stderr, "Error: -strand must be -1, 1, or 2\n");
        return 1;
    }

    if (config.mode < 1 || config.mode > 3) {
        std::fprintf(stderr, "Error: -mode must be 1, 2, or 3\n");
        return 1;
    }

    // Mode x stage consistency: an option that only takes effect in a later
    // stage is an error when an earlier mode is selected.  Stage 2 options
    // need mode >= 2; Stage 3 options need mode 3.
    {
        static const char* const kStage2Opts[] = {
            "-stage2_min_score", "-stage2_max_gap", "-stage2_max_lookback",
            "-stage2_min_nhit_diag", "-stage2_max_nhit_per_subject",
            "-stage2_max_nhit_per_subject_mode", "-stage2_max_nhit_in_total",
        };
        static const char* const kStage3Opts[] = {
            "-stage3_traceback", "-stage3_gapopen", "-stage3_gapext",
            "-stage3_min_ppositive", "-stage3_min_npositive",
            "-stage3_score_matrix", "-context_extend", "-db",
            "-stage3_max_nhit_per_subject", "-stage3_max_nhit_per_subject_mode",
            "-stage3_max_nhit_in_total",
        };
        if (config.mode == 1) {
            for (const char* opt : kStage2Opts) {
                if (cli.has(opt)) {
                    std::fprintf(stderr,
                        "Error: %s requires -mode 2 or higher\n", opt);
                    return 1;
                }
            }
            for (const char* opt : kStage3Opts) {
                if (cli.has(opt)) {
                    std::fprintf(stderr, "Error: %s requires -mode 3\n", opt);
                    return 1;
                }
            }
        } else if (config.mode == 2) {
            for (const char* opt : kStage3Opts) {
                if (cli.has(opt)) {
                    std::fprintf(stderr, "Error: %s requires -mode 3\n", opt);
                    return 1;
                }
            }
        }
    }

    // BLAST DB path for mode 3 (default: same as index prefix)
    std::string db_path = cli.get_string("-db", ix_prefix);

    // Stage 3 config
    Stage3Config stage3_config;
    stage3_config.gapopen = cli.get_int("-stage3_gapopen", 10);
    stage3_config.gapext = cli.get_int("-stage3_gapext", 1);
    stage3_config.traceback = (cli.get_int("-stage3_traceback", 0) != 0);
    stage3_config.min_ppositive = cli.get_double("-stage3_min_ppositive", 0.0);
    stage3_config.min_npositive = static_cast<uint32_t>(cli.get_int("-stage3_min_npositive", 0));
    stage3_config.max_nhit_per_subject =
        static_cast<uint32_t>(cli.get_int("-stage3_max_nhit_per_subject", 1));
    if (cli.has("-stage3_max_nhit_per_subject_mode")) {
        int m = cli.get_int("-stage3_max_nhit_per_subject_mode", 3);
        if (m < 1 || m > 4) {
            std::fprintf(stderr,
                "Error: -stage3_max_nhit_per_subject_mode must be 1, 2, 3, or 4\n");
            return 1;
        }
        stage3_config.max_nhit_per_subject_mode = static_cast<uint8_t>(m);
    }
    stage3_config.max_nhit_in_total =
        static_cast<uint32_t>(cli.get_int("-stage3_max_nhit_in_total", 0));
    {
        std::string err;
        if (!parse_score_matrix(cli, stage3_config.score_matrix, err)) {
            std::fprintf(stderr, "%s\n", err.c_str());
            return 1;
        }
    }
    // Parse -context_extend: integer (bases) or decimal (query length multiplier)
    ContextParam ctx_param;
    {
        std::string err;
        if (!parse_context(cli.get_string("-context_extend", "2.0"), ctx_param, err)) {
            std::fprintf(stderr, "%s\n", err.c_str());
            return 1;
        }
    }

    // Parse -stage1_min_score as double to support fractional values
    {
        double min_s1 = cli.get_double("-stage1_min_score", 0.5);
        if (min_s1 > 0 && min_s1 < 1.0) {
            config.min_stage1_score_frac = min_s1;
            // Leave stage1.min_stage1_score at default (will be overridden per-query)
        } else {
            config.stage1.min_stage1_score = static_cast<uint32_t>(min_s1);
        }
    }

    config.accept_qdegen = static_cast<uint8_t>(cli.get_int("-accept_qdegen", 1));

    // -min_query_length default 64.  The integrity check
    // against the index's min_seq_length runs after the first .kix is
    // opened (search.cpp keeps a single in-memory list of probes).
    {
        int v = cli.get_int("-min_query_length", 64);
        if (v <= 0) {
            std::fprintf(stderr,
                "Error: -min_query_length must be a positive integer (got %d)\n", v);
            return 1;
        }
        config.min_query_length = static_cast<uint32_t>(v);
    }

    {
        std::string err;
        if (!parse_max_degen_expand(cli, 16, config.max_degen_expand, err)) {
            std::fprintf(stderr, "%s\n", err.c_str());
            return 1;
        }
    }

    uint8_t spaced_t = 0;
    {
        std::string err;
        if (!parse_spaced_seed_t(cli, spaced_t, err)) {
            std::fprintf(stderr, "%s\n", err.c_str());
            return 1;
        }
    }

    TemplateType spaced_type = TemplateType::kBoth;
    {
        std::string err;
        if (!parse_template_type_cli(cli, TemplateType::kBoth,
                                     /*allow_contiguous=*/false,
                                     spaced_type, err)) {
            std::fprintf(stderr, "%s\n", err.c_str());
            return 1;
        }
    }

    // When -k is given, validate the (k, t) pair up front; when -k is omitted
    // the pair is re-validated against the resolved variant's k below.
    if (spaced_t > 0 && filter_k > 0) {
        if (!validate_spaced_seed(filter_k, spaced_t)) {
            std::fprintf(stderr, "Error: -t %d is not valid for -k %d\n", spaced_t, filter_k);
            return 1;
        }
    }

    config.t = spaced_t;

    // Output format
    OutputFormat outfmt;
    {
        std::string err;
        if (!parse_output_format(cli.get_string("-output_format", "tsv"), outfmt, err)) {
            std::fprintf(stderr, "%s\n", err.c_str());
            return 1;
        }
    }

    // Validate output format compatibility
    {
        std::string err;
        if (!validate_output_format(outfmt, config.mode, stage3_config.traceback,
                                    output_path, err)) {
            std::fprintf(stderr, "%s\n", err.c_str());
            return 1;
        }
    }

    // Validate compression level (parsed against the codec selected by the
    // output path's trailing suffix).  When -compression_level is absent we
    // pass the kCompressionLevelDefault sentinel through unchanged.
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
    // mode 3 requires BLAST DB
    if (config.mode == 3) {
        auto vol_paths = BlastDbReader::find_volume_paths(db_path);
        if (vol_paths.empty()) {
            std::fprintf(stderr, "Error: -mode 3 requires BLAST DB but none found at '%s'\n",
                         db_path.c_str());
            return 1;
        }
    }

    // Discover every variant under the prefix, then narrow to a single one
    // with the shared variant filter.  The template_type dimension is resolved
    // per the agreed semantics:
    //   -t 0                       -> single contiguous
    //   -t>0 -template_type cod/opt -> single coding / optimal
    //   -t>0 -template_type both    -> coding+optimal pair (both required)
    //   -t>0 (no -template_type)    -> pair if both exist, else the present side
    TemplateResolveMode resolve_mode;
    if (spaced_t == 0) {
        resolve_mode = TemplateResolveMode::kSingle;            // contiguous
    } else if (!cli.has("-template_type")) {
        resolve_mode = TemplateResolveMode::kBothOrSingle;      // wildcard
    } else if (spaced_type == TemplateType::kBoth) {
        resolve_mode = TemplateResolveMode::kBothRequired;
    } else {
        resolve_mode = TemplateResolveMode::kSingle;            // coding / optimal
    }

    VariantFilter vf;
    vf.has_k = (filter_k > 0); vf.k = filter_k;
    vf.t_is_wildcard = false;  vf.t = spaced_t;   // unset -t defaults to 0
    if (cli.has("-template_type")) {
        vf.has_template_type = true;
        vf.template_type = static_cast<uint8_t>(spaced_type);
    }
    if (cli.has("-min_seq_length")) {
        vf.has_min_seq_length = true;
        vf.min_seq_length = static_cast<uint32_t>(cli.get_int("-min_seq_length", 0));
    }
    if (cli.has("-min_length_split")) {
        vf.has_min_length_split = true;
        vf.min_length_split = static_cast<uint32_t>(cli.get_int("-min_length_split", 0));
    }
    if (cli.has("-overlap_length")) {
        vf.has_overlap_length = true;
        vf.overlap_length = static_cast<uint32_t>(cli.get_int("-overlap_length", 0));
    }
    if (cli.has("-max_freq_build")) {
        vf.has_max_freq_build = true;
        vf.max_freq_build = std::strtoull(
            cli.get_string("-max_freq_build").c_str(), nullptr, 10);
    }
    // -max_degen_expand is a query parameter (config.max_degen_expand) and the
    // variant-selection tie-break, never a filter.

    // Pre-narrow by k at discovery (cheap, and avoids stat()ing other-k
    // volumes); variant_matches() below applies the remaining filters.
    auto all_vols = discover_volumes(ix_prefix, filter_k);
    std::vector<DiscoveredVolume> matched;
    std::vector<VariantFields> matched_fields;
    for (const auto& v : all_vols) {
        VariantFields f = variant_fields_of(v);
        if (variant_matches(vf, f)) {
            matched.push_back(v);
            matched_fields.push_back(f);
        }
    }

    std::vector<TemplateContext> ctxs;
    ctxs.reserve(2);  // pin addresses; khx_ptr depends on stable storage
    bool is_both_mode = false;
    {
        std::string err;
        auto sel = resolve_variant_indices(matched_fields, resolve_mode,
                                           config.max_degen_expand, err);
        if (sel.empty()) {
            std::fprintf(stderr, "Error: %s (prefix %s)\n",
                         err.c_str(), ix_prefix.c_str());
            return 1;
        }
        // Infer single-vs-pair from the resolved template_types.
        bool has_cod = false, has_opt = false;
        for (size_t i : sel) {
            if (matched[i].template_type == 1) has_cod = true;
            else if (matched[i].template_type == 2) has_opt = true;
        }
        is_both_mode = has_cod && has_opt;

        if (is_both_mode) {
            ctxs.emplace_back();
            ctxs.back().type = TemplateType::kCoding;
            ctxs.emplace_back();
            ctxs.back().type = TemplateType::kOptimal;
            for (size_t i : sel) {
                if (matched[i].template_type == 1)
                    ctxs[0].vol_files.push_back(matched[i]);
                else
                    ctxs[1].vol_files.push_back(matched[i]);
            }
            std::sort(ctxs[0].vol_files.begin(), ctxs[0].vol_files.end(),
                      [](const DiscoveredVolume& a, const DiscoveredVolume& b) {
                          return a.volume_index < b.volume_index;
                      });
            std::sort(ctxs[1].vol_files.begin(), ctxs[1].vol_files.end(),
                      [](const DiscoveredVolume& a, const DiscoveredVolume& b) {
                          return a.volume_index < b.volume_index;
                      });
            if (ctxs[0].vol_files.size() != ctxs[1].vol_files.size()) {
                std::fprintf(stderr,
                    "Error: coding and optimal volume counts differ (%zu vs %zu)\n",
                    ctxs[0].vol_files.size(), ctxs[1].vol_files.size());
                return 1;
            }
        } else {
            // Single variant: every selected volume shares one template_type,
            // so take it from the resolved variant (contiguous => kContiguous).
            ctxs.emplace_back();
            ctxs.back().type =
                static_cast<TemplateType>(matched[sel[0]].template_type);
            for (size_t i : sel) ctxs.back().vol_files.push_back(matched[i]);
        }
    }

    // ctxs[0] is the canonical side (coding for both-mode, the only side otherwise).
    const auto& vol_files = ctxs[0].vol_files;

    // The selected variant is unique in k, so read it straight off.
    const int k = vol_files.front().k;

    if (spaced_t > 0 && !validate_spaced_seed(k, spaced_t)) {
        std::fprintf(stderr, "Error: -t %d is not valid for k=%d\n", spaced_t, k);
        return 1;
    }

    logger.info("Found %zu volume(s), k=%d, threads=%d", vol_files.size(), k, nthread);

    // Read query FASTA or primer FASTA
    std::vector<FastaRecord> queries;
    std::vector<PrimerPair> primer_pairs;

    // Primer mode parameters
    uint32_t insert_length = 0;
    double stage1_primer_score = 0.5;
    int stage2_primer_score_add = 1;

    if (has_primer) {
        insert_length = static_cast<uint32_t>(cli.get_int("-insert_length", 0));
        stage1_primer_score = cli.get_double("-stage1_primer_score", 0.5);
        stage2_primer_score_add = cli.get_int("-stage2_primer_score_add", 1);

        if (stage1_primer_score > 1.0 && stage1_primer_score < 2.0) {
            std::fprintf(stderr, "Error: -stage1_primer_score must be 0<v<=1 (fraction) or v>=2 (absolute)\n");
            return 1;
        }

        auto primer_records = read_fasta(primer_path);
        if (primer_records.empty()) {
            std::fprintf(stderr, "Error: no primer sequences found\n");
            return 1;
        }

        // Park the raw primer records here; they are reparsed into real
        // queries below once k and seed_masks are resolved.
        queries.reserve(primer_records.size());
        for (auto& r : primer_records) {
            queries.push_back(std::move(r));
        }
        logger.info("Read %zu primer sequence(s)", queries.size());
    } else {
        queries = read_fasta(query_path);
        if (queries.empty()) {
            std::fprintf(stderr, "Error: no query sequences found\n");
            return 1;
        }
        logger.info("Read %zu query sequence(s)", queries.size());
    }

    // skip_reason / skip_detail per query (filled in by preprocess_query below).
    std::vector<uint8_t> query_skip_reason(queries.size(), 0);
    std::vector<std::string> query_skip_detail(queries.size());
    bool has_skipped = false;

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

    // Helper lambda: validate every kix/kpx file and capture posting / full
    // sizes into VolumeMeta.  Readers are opened once for header validation
    // and immediately closed; the search orchestrator re-opens kix/kpx
    // per-batch with explicit MADV_WILLNEED / MADV_RANDOM (.ksx stays open
    // for accession / seq_length lookups and OidFilter construction).
    auto open_volumes = [&](const std::vector<DiscoveredVolume>& vfiles,
                            std::vector<VolumeData>& vdata,
                            bool need_kpx) -> bool {
        vdata.resize(vfiles.size());
        for (size_t vi = 0; vi < vfiles.size(); vi++) {
            const auto& vf = vfiles[vi];
            vdata[vi].files = vf;
            vdata[vi].volume_index = vf.volume_index;

            // Validate kix and capture num_sequences.
            {
                KixReader kix_probe;
                if (!kix_probe.open(vf.kix_path)) {
                    std::fprintf(stderr, "Error: cannot open %s\n",
                                 vf.kix_path.c_str());
                    return false;
                }
                vdata[vi].num_sequences    = kix_probe.num_sequences();
                // -min_query_length must be >= the index's min_seq_length.
                const uint32_t idx_min = vf.min_seq_length;
                if (config.min_query_length < idx_min) {
                    std::fprintf(stderr,
                        "Error: -min_query_length=%u is smaller than the "
                        "index's min_seq_length=%u (%s). "
                        "Specify -min_query_length=%u or larger.\n",
                        static_cast<unsigned>(config.min_query_length),
                        static_cast<unsigned>(idx_min),
                        vf.kix_path.c_str(),
                        static_cast<unsigned>(idx_min));
                    kix_probe.close();
                    return false;
                }
                kix_probe.close();
            }

            if (need_kpx) {
                if (!vf.has_kpx) {
                    std::fprintf(stderr,
                        "Error: mode %d requires .kpx files, but %s was not found.\n"
                        "This index was built with -mode 1 (stage 1 only).\n",
                        config.mode, vf.kpx_path.c_str());
                    return false;
                }
                KpxReader kpx_probe;
                if (!kpx_probe.open(vf.kpx_path)) {
                    std::fprintf(stderr, "Error: cannot open %s\n",
                                 vf.kpx_path.c_str());
                    return false;
                }
                kpx_probe.close();
            }

            if (!vdata[vi].ksx.open(vf.ksx_path)) {
                std::fprintf(stderr, "Error: cannot open %s\n", vf.ksx_path.c_str());
                return false;
            }
            if (filter_mode != OidFilterMode::kNone) {
                vdata[vi].filter.build(seqidlist, vdata[vi].ksx, filter_mode);
            }
        }
        return true;
    };

    const bool need_kpx = (config.mode != 1);

    // Pre-open / validate volumes for every context.
    for (auto& ctx : ctxs) {
        if (!open_volumes(ctx.vol_files, ctx.volumes, need_kpx)) return 1;
    }

    // Drive max_query_length from the index's overlap_length; 0
    // disables the upper-bound check (no fragment splitting).
    if (!ctxs.empty() && !ctxs[0].vol_files.empty()) {
        config.max_query_length = ctxs[0].vol_files[0].overlap_length;
    }

    // Open shared .khx for every context (non-fatal if missing).
    {
        auto parts = parse_index_prefix(ix_prefix);
        for (auto& ctx : ctxs) {
            if (ctx.vol_files.empty()) continue;
            const auto& v0 = ctx.vol_files[0];
            ctx.khx.open(khx_path_for(parts.parent_dir, parts.db, k,
                                       spaced_t, static_cast<uint8_t>(ctx.type),
                                       v0.min_seq_length, v0.min_length_split,
                                       v0.overlap_length, v0.max_freq_build,
                                       v0.max_degen_expand));
        }
    }

    // Persistent madvise covers only the small per-volume metadata that
    // every search re-walks (.khx + .ksx).  The .kix / .kpx mappings live
    // for one group at a time inside the search orchestrator, so they are
    // not pinned here.  Stage 1 / Stage 2A no longer spend a memory budget;
    // whatever budget remains after .khx/.ksx becomes `posting_budget`,
    // which bounds the Stage 3 posting heap batching.
    uint64_t posting_budget = 0;
    {
        uint64_t budget = memory_limit;
        auto try_willneed = [&budget](auto& reader) {
            uint64_t sz = reader.willneed_size();
            bool fits = (sz > 0 && budget >= sz);
            reader.apply_madvise(fits);
            if (fits) budget -= sz;
        };
        for (auto& ctx : ctxs) try_willneed(ctx.khx);
        for (auto& vd : ctxs[0].volumes) try_willneed(vd.ksx);
        posting_budget = budget;
        logger.info("madvise budget: %s used / %s total (posting_budget=%s)",
                    format_size(memory_limit - budget).c_str(),
                    format_size(memory_limit).c_str(),
                    format_size(posting_budget).c_str());
    }

    // Per-context derived state: KHX pointer, seed masks.
    for (auto& ctx : ctxs) {
        ctx.khx_ptr = ctx.khx.is_open() ? &ctx.khx : nullptr;
        if (spaced_t > 0) ctx.seed_masks = get_seed_masks(k, spaced_t, ctx.type);
    }

    // Primer mode: parse pairs and generate query sequences (needs k and seed_masks)
    if (has_primer) {
        auto primer_records = std::move(queries);
        queries.clear();

        PrimerConfig pcfg;
        pcfg.insert_length = insert_length;
        pcfg.k = k;
        pcfg.t = spaced_t;
        // For primer mode: use ctxs[0]'s masks (coding side in both-mode) for
        // position counting.
        pcfg.masks = spaced_t > 0 ? &ctxs[0].seed_masks : nullptr;

        std::string primer_err = parse_primer_pairs(primer_records, pcfg, primer_pairs);
        if (!primer_err.empty()) {
            std::fprintf(stderr, "%s\n", primer_err.c_str());
            return 1;
        }

        // Convert primer pairs to queries and resolve thresholds
        uint32_t min_stage2_score = UINT32_MAX;
        for (const auto& pp : primer_pairs) {
            FastaRecord qr;
            qr.id = pp.query_id;
            qr.sequence = pp.query_seq;
            queries.push_back(std::move(qr));

            uint32_t total_pos = pp.fwd_kmer_positions + pp.rev_kmer_positions;

            // Validate stage1_primer_score (absolute mode)
            if (stage1_primer_score >= 2.0) {
                if (static_cast<uint32_t>(stage1_primer_score) > total_pos) {
                    std::fprintf(stderr,
                        "Error: -stage1_primer_score (%u) exceeds total k-mer positions (%u) for pair %s\n",
                        static_cast<unsigned>(stage1_primer_score), total_pos, pp.query_id.c_str());
                    return 1;
                }
            }

            // Compute stage2 threshold: max(fwd_pos, rev_pos) + add
            uint32_t s2_score = std::max(pp.fwd_kmer_positions, pp.rev_kmer_positions)
                                + static_cast<uint32_t>(stage2_primer_score_add);
            if (s2_score > total_pos) {
                std::fprintf(stderr,
                    "Error: stage2 threshold (%u) exceeds total k-mer positions (%u) for pair %s\n",
                    s2_score, total_pos, pp.query_id.c_str());
                return 1;
            }
            min_stage2_score = std::min(min_stage2_score, s2_score);
        }

        // Set stage1 threshold
        if (stage1_primer_score > 0 && stage1_primer_score <= 1.0) {
            config.min_stage1_score_frac = stage1_primer_score;
        } else if (stage1_primer_score >= 2.0) {
            config.stage1.min_stage1_score = static_cast<uint32_t>(stage1_primer_score);
        }

        // Set stage2 threshold (use minimum across all pairs)
        config.stage2.min_score = min_stage2_score;

        // Set stage2 max_gap to insert_length (unless explicitly overridden)
        if (!cli.has("-stage2_max_gap")) {
            config.stage2.max_gap = insert_length;
        }

        // Reset skip tracking for the new queries
        query_skip_reason.assign(queries.size(), 0);
        query_skip_detail.assign(queries.size(), std::string());
        has_skipped = false;

        logger.info("Generated %zu primer pair queries", primer_pairs.size());
    }

    // Preprocess queries in parallel.  Slots are pre-allocated so each
    // task writes to ctx.pp{16,32}[qi] independently; query_pp_idx
    // collapses to the identity (qi -> qi).  Warnings and the
    // skip-reason bookkeeping run in a sequential post-pass for
    // deterministic ordering.
    const bool use_uint16 = (kmer_type_for(k, spaced_t) == 0);
    std::vector<size_t> query_pp_idx(queries.size(), SIZE_MAX);
    for (size_t qi = 0; qi < queries.size(); ++qi) query_pp_idx[qi] = qi;

    for (auto& ctx : ctxs) {
        if (use_uint16) ctx.pp16.resize(queries.size());
        else            ctx.pp32.resize(queries.size());
    }

    std::vector<uint8_t> any_multi_degen_per_q(queries.size(), 0);

    {
        tbb::task_arena arena_pp(nthread);
        arena_pp.execute([&] {
            tbb::parallel_for(
                tbb::blocked_range<size_t>(0, queries.size()),
                [&](const tbb::blocked_range<size_t>& range) {
                    for (size_t qi = range.begin(); qi != range.end(); ++qi) {
                        if (use_uint16) {
                            std::vector<PreprocessTemplateBinding<uint16_t>> b(ctxs.size());
                            for (size_t i = 0; i < ctxs.size(); ++i) {
                                b[i].slot  = &ctxs[i].pp16[qi].qdata;
                                b[i].khx   = ctxs[i].khx_ptr;
                                b[i].masks = &ctxs[i].seed_masks;
                            }
                            auto out = run_preprocess_one_query<uint16_t>(
                                queries[qi].sequence, k, spaced_t, config,
                                b.data(), b.size());
                            if (out.skip_reason != 0) {
                                query_skip_reason[qi] = out.skip_reason;
                                query_skip_detail[qi] = std::move(out.skip_detail);
                            }
                            any_multi_degen_per_q[qi] = out.multi_degen ? 1 : 0;
                        } else {
                            std::vector<PreprocessTemplateBinding<uint32_t>> b(ctxs.size());
                            for (size_t i = 0; i < ctxs.size(); ++i) {
                                b[i].slot  = &ctxs[i].pp32[qi].qdata;
                                b[i].khx   = ctxs[i].khx_ptr;
                                b[i].masks = &ctxs[i].seed_masks;
                            }
                            auto out = run_preprocess_one_query<uint32_t>(
                                queries[qi].sequence, k, spaced_t, config,
                                b.data(), b.size());
                            if (out.skip_reason != 0) {
                                query_skip_reason[qi] = out.skip_reason;
                                query_skip_detail[qi] = std::move(out.skip_detail);
                            }
                            any_multi_degen_per_q[qi] = out.multi_degen ? 1 : 0;
                        }
                    }
                });
        });
    }

    uint32_t skipped_too_short = 0;
    uint32_t skipped_too_long  = 0;
    for (size_t qi = 0; qi < queries.size(); ++qi) {
        if (query_skip_reason[qi] != 0) {
            has_skipped = true;
            if (query_skip_reason[qi] == kSkipQueryTooShort) ++skipped_too_short;
            else if (query_skip_reason[qi] == kSkipQueryTooLong) ++skipped_too_long;
            std::fprintf(stderr, "Warning: query '%s' skipped: %s (%s)\n",
                         queries[qi].id.c_str(),
                         skip_reason_str(query_skip_reason[qi]),
                         query_skip_detail[qi].c_str());
        }
        if (any_multi_degen_per_q[qi]) {
            std::fprintf(stderr,
                "Warning: query '%s' contains k-mers exceeding max_degen_expand=%u; "
                "those k-mers are ignored and not used in the search\n",
                queries[qi].id.c_str(),
                static_cast<unsigned>(config.max_degen_expand));
        }
    }
    if (skipped_too_short > 0) {
        logger.info("Skipped %u queries shorter than min_query_length=%u",
                    skipped_too_short,
                    static_cast<unsigned>(config.min_query_length));
    }
    if (skipped_too_long > 0) {
        logger.info("Skipped %u queries longer than overlap_length=%u",
                    skipped_too_long,
                    static_cast<unsigned>(config.max_query_length));
    }

    // Thread-local Stage1Buffer sizing inputs.  num_sequences was captured
    // at validation time; width is chosen from actual preprocessed k-mer
    // counts.  The orchestrator constructs the TLS buffer pool internally.
    uint32_t max_num_seqs = 0;
    for (const auto& vd : ctxs[0].volumes)
        max_num_seqs = std::max(max_num_seqs, vd.num_sequences);

    uint32_t max_kmer_positions = 0;
    for (const auto& ctx : ctxs) {
        max_kmer_positions = use_uint16
            ? accumulate_max_kmer_positions(max_kmer_positions, ctx.pp16)
            : accumulate_max_kmer_positions(max_kmer_positions, ctx.pp32);
    }
    Stage1Width width = select_width(max_kmer_positions, max_kmer_positions);

    size_t num_volumes = ctxs[0].volumes.size();

    // Build VolumeMeta + ksx / OidFilter inputs and dispatch through the
    // common search orchestrator.  ctxs[*].volumes outlive the call.
    auto run_orchestrated = [&](auto kmer_int_tag) {
        using KmerInt = decltype(kmer_int_tag);

        RunSearchInputs<KmerInt> in;
        in.both_mode = is_both_mode;
        in.k = k;
        in.nthread = nthread;
        in.config = config;
        in.logger = &logger;
        in.max_num_seqs = max_num_seqs;
        in.width = width;

        in.volumes_cod.resize(num_volumes);
        if (is_both_mode) in.volumes_opt.resize(num_volumes);
        in.ksx_per_volume.resize(num_volumes);
        in.oid_filters.resize(num_volumes);

        for (size_t vi = 0; vi < num_volumes; ++vi) {
            const auto& v0 = ctxs[0].volumes[vi];
            in.volumes_cod[vi].files             = v0.files;
            in.volumes_cod[vi].volume_index      = v0.volume_index;
            in.volumes_cod[vi].num_sequences     = v0.num_sequences;
            if (is_both_mode) {
                const auto& v1 = ctxs[1].volumes[vi];
                in.volumes_opt[vi].files             = v1.files;
                in.volumes_opt[vi].volume_index      = v1.volume_index;
                in.volumes_opt[vi].num_sequences     = v1.num_sequences;
            }
            in.ksx_per_volume[vi] = &ctxs[0].volumes[vi].ksx;
            in.oid_filters[vi]    = std::move(ctxs[0].volumes[vi].filter);
        }

        std::vector<QueryBundle<KmerInt>> query_bundles(queries.size());
        for (size_t qi = 0; qi < queries.size(); ++qi) {
            query_bundles[qi].query_id = &queries[qi].id;
            if (query_skip_reason[qi] != 0) continue;
            size_t pp_idx = query_pp_idx[qi];
            if constexpr (std::is_same_v<KmerInt, uint16_t>) {
                query_bundles[qi].qdata_primary = &ctxs[0].pp16[pp_idx].qdata;
                if (is_both_mode)
                    query_bundles[qi].qdata_secondary = &ctxs[1].pp16[pp_idx].qdata;
            } else {
                query_bundles[qi].qdata_primary = &ctxs[0].pp32[pp_idx].qdata;
                if (is_both_mode)
                    query_bundles[qi].qdata_secondary = &ctxs[1].pp32[pp_idx].qdata;
            }
        }

        in.queries = &query_bundles;
        in.query_skip_reason = &query_skip_reason;

        return run_search<KmerInt>(in);
    };

    std::vector<OrchestratorHit> orch_hits;
    if (use_uint16) orch_hits = run_orchestrated(uint16_t{});
    else            orch_hits = run_orchestrated(uint32_t{});

    // ----------------------------------------------------------------
    // Mode 1 parallel TSV / JSON path.
    //
    // Mode 1's hit dump is flat (no Stage 3, no SAM) and the
    // OrchestratorHit -> OutputHit conversion + serial std::ostream
    // formatting is the dominant single-thread tail of the run on
    // large databases.  When the output format is TSV or JSON we
    // route OrchestratorHit straight into the parallel formatter
    // and skip the OutputHit pipeline entirely.  Mode 1 SAM/BAM is
    // rejected at validate_output_format(), and Mode 2/3 still need
    // the OutputHit path for chain / Stage 3 / SAM fields.
    // ----------------------------------------------------------------
    const bool use_parallel_mode1 = (config.mode == 1) &&
        (outfmt == OutputFormat::kTsv || outfmt == OutputFormat::kJson);

    if (use_parallel_mode1) {
        // Sort by query_idx (JSON groups by query, TSV preserves input order).
        // The per-query result cap is the Stage 1 in-total limit, already
        // applied inside run_search.
        tbb::parallel_sort(orch_hits.begin(), orch_hits.end(),
            [](const OrchestratorHit& a, const OrchestratorHit& b) {
                return a.query_idx < b.query_idx;
            });

        // Open compressed output sink.
        std::string err;
        auto out = open_output_compressed(output_path, compression_level, err);
        if (!out) {
            std::fprintf(stderr, "%s\n", err.c_str());
            return 1;
        }

        Mode1ParallelInputs in;
        in.hits        = &orch_hits;
        in.queries     = &queries;
        in.skip_reason = &query_skip_reason;
        in.skip_detail = &query_skip_detail;
        in.nthread     = nthread;
        in.ksx_per_volume.resize(num_volumes);
        for (size_t vi = 0; vi < num_volumes; ++vi)
            in.ksx_per_volume[vi] = &ctxs[0].volumes[vi].ksx;

        logger.info("Writing %zu hit(s)...", orch_hits.size());
        if (outfmt == OutputFormat::kTsv)
            write_results_tsv_mode1_parallel(*out.stream, in);
        else
            write_results_json_mode1_parallel(*out.stream, in);
        logger.info("Done. %zu hit(s) reported.", orch_hits.size());
        return has_skipped ? 2 : 0;
    }

    // Convert OrchestratorHit -> OutputHit at the boundary.
    //
    // Stage 2 chains live in fragment-relative coordinates.
    // Re-map them to parent-relative coordinates for output: sseqid
    // becomes the parent accession, slen becomes the parent length, and
    // sstart/send shift by (fragment_start - 1).
    std::vector<OutputHit> all_hits;
    all_hits.reserve(orch_hits.size());
    for (const auto& oh_in : orch_hits) {
        const auto& cr = oh_in.cr;
        const auto& ksx_primary = ctxs[0].volumes[oh_in.volume_idx].ksx;
        const uint32_t parent_idx = ksx_primary.parent_index(cr.seq_id);
        const uint32_t frag_start = ksx_primary.fragment_start(cr.seq_id);
        const uint32_t shift      = frag_start - 1u;
        OutputHit oh;
        oh.qseqid = queries[oh_in.query_idx].id;
        oh.sseqid = std::string(ksx_primary.parent_accession(parent_idx));
        oh.sstrand = cr.is_reverse ? '-' : '+';
        oh.qstart = cr.q_start;
        oh.qend = cr.q_end;
        oh.sstart = cr.s_start + shift;
        oh.send   = cr.s_end   + shift;
        oh.chainscore = cr.chainscore;
        oh.coverscore = cr.stage1_score;
        oh.volume = oh_in.volume_index;
        // Stage 3 fetches the subject sequence via BlastDbReader keyed by
        // BLAST DB volume-local OID, so the parent's BLAST OID is what
        // belongs in oh.oid (not the internal fragment seq_id).
        oh.oid = ksx_primary.blast_oid(parent_idx);
        oh.qlen = static_cast<uint32_t>(queries[oh_in.query_idx].sequence.size());
        oh.slen = ksx_primary.parent_length(parent_idx);
        all_hits.push_back(std::move(oh));
    }

    // Build skip-marker rows for queries that were not searched. Held aside
    // so Stage 3 only sees real hits.
    std::vector<OutputHit> skip_markers;
    for (size_t qi = 0; qi < queries.size(); qi++) {
        if (query_skip_reason[qi] == 0) continue;
        OutputHit oh;
        oh.qseqid = queries[qi].id;
        oh.qlen = static_cast<uint32_t>(queries[qi].sequence.size());
        oh.skip_reason = query_skip_reason[qi];
        oh.skip_detail = query_skip_detail[qi];
        oh.sstrand = '*';
        skip_markers.push_back(std::move(oh));
    }

    // Stage 3 alignment (mode 3 only)
    if (config.mode == 3) {
        logger.info("Running Stage 3 alignment on %zu hits...", all_hits.size());
        stage3_config.posting_budget = posting_budget;
        all_hits = run_stage3(all_hits, queries, db_path, stage3_config,
                              ctx_param.is_ratio, ctx_param.ratio, ctx_param.abs,
                              logger);
        logger.info("Stage 3 complete: %zu hits after filtering.", all_hits.size());
        // Stage 3 dedup over parent-relative (send, alnscore).
        // Re-runs after alignment because the alignment can extend / clamp
        // sstart freely, leaving duplicates that survived Stage 2 dedup with
        // distinct (sstart, chainscore) but identical (send, alnscore).
        const size_t before = all_hits.size();
        dedup_stage3_output_hits(all_hits);
        if (before != all_hits.size()) {
            logger.info("Stage 3 dedup: %zu hit(s) -> %zu after dedup",
                        before, all_hits.size());
        }

        // Stage 3 per-subject (N) and in-total (L) caps, by alnscore.  Applied
        // after dedup so duplicates do not count against the limits.
        select_parent_topn_output(all_hits,
                                  stage3_config.max_nhit_per_subject,
                                  stage3_config.max_nhit_per_subject_mode);
        apply_in_total_output(all_hits, stage3_config.max_nhit_in_total);
    }

    // Re-attach skip markers after Stage 3.
    for (auto& m : skip_markers) all_hits.push_back(std::move(m));

    // Write output
    if (!write_all_results(output_path, all_hits, outfmt,
                           config.mode,
                           stage3_config.traceback, compression_level)) {
        return 1;
    }

    logger.info("Done. %zu hit(s) reported.", all_hits.size());
    return has_skipped ? 2 : 0;
}
