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
#include "search/volume_searcher.hpp"
#include "search/parallel_search.hpp"
#include "search/preprocess_runner.hpp"
#include "search/query_preprocessor.hpp"
#include "search/stage3_alignment.hpp"
#include "search/tier_selection.hpp"
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
#include <set>
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
        "  -memory_limit <size>     madvise WILLNEED budget (default: half of RAM)\n"
        "                           Accepts K, M, G suffixes\n"
        "  -mode <1|2|3>            1=Stage1, 2=Stage1+2, 3=Stage1+2+3 (default: 2)\n"
        "  -db <path>               BLAST DB path for mode 3 (default: same as -ix)\n"
        "  -stage1_score <1|2>      1=coverscore, 2=matchscore (default: 1)\n"
        "  -stage2_min_score <int>  Minimum chain score (default: 0 = adaptive)\n"
        "                           0 = use resolved Stage 1 threshold\n"
        "  -stage2_max_gap <int>    Chaining diagonal gap tolerance (default: 100)\n"
        "  -stage2_max_lookback <int>  Chaining DP lookback window (default: 64, 0=unlimited)\n"
        "  -stage2_max_nhit_per_subject <int>  Max chains per subject (default: 1, 0=unlimited)\n"
        "  -stage2_min_nhit_diag <int>  Diagonal filter min hits (default: 1)\n"
        "  -stage1_topn <int>       Stage 1 candidate limit, 0=unlimited (default: 0)\n"
        "  -stage1_min_score <num>  Stage 1 minimum score; integer or 0<P<1 fraction (default: 0.5)\n"
        "  -nresult <int>           Max results per query, 0=unlimited (default: 0)\n"
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
        "  -stage3_score_matrix <name>  Score matrix: degmatch, dnafull, nuc44 (default: degmatch)\n"
        "  -stage3_nthread_fetch <int>  Threads for BLAST DB fetch in mode 3 (default: min(8, nthread))\n"
        "  -max_degen_expand <int>  Max degenerate expansion per k-mer (default: 16, max: 256, 0/1: disable)\n"
        "  -t <int>                 Template length for spaced seeds (0=contiguous, 13/15/18 for k=8-9, 16/18/21 for k=11-12; default: 0)\n"
        "  -template_type <string>  Template type: coding, optimal, both (default: both)\n"
        "  -output_format <tsv|json|sam|bam>  Output format (default: tsv)\n"
        "  -compression_level <int>    Output compression level (default per codec: gzip=6, bzip2=9, xz=6, zstd=3)\n"
        "  -v, --verbose            Verbose logging\n",
        prog);
}

// Pre-opened volume data shared across threads (read-only after init).
struct VolumeData {
    KixReader kix;
    KpxReader kpx;
    KsxReader ksx;
    OidFilter filter;
    uint16_t volume_index;
};

// Per-template-side data. For "both" template_type the search holds two
// contexts (coding + optimal); otherwise a single context is used. Each
// context owns its own volume readers, KHX bitset, seed masks, and
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
    config.stage1.stage1_topn = static_cast<uint32_t>(cli.get_int("-stage1_topn", 0));
    config.stage1.stage1_score_type = static_cast<uint8_t>(cli.get_int("-stage1_score", 1));
    config.stage2.max_gap = static_cast<uint32_t>(cli.get_int("-stage2_max_gap", 100));
    config.stage2.chain_max_lookback = static_cast<uint32_t>(cli.get_int("-stage2_max_lookback", 64));
    config.stage2.max_nhit_per_subject = static_cast<uint32_t>(cli.get_int("-stage2_max_nhit_per_subject", 1));
    config.stage2.min_nhit_diag = static_cast<uint32_t>(cli.get_int("-stage2_min_nhit_diag", 1));
    config.stage2.min_score = static_cast<uint32_t>(cli.get_int("-stage2_min_score", 0));
    config.nresult = static_cast<uint32_t>(cli.get_int("-nresult", 0));
    config.mode = static_cast<uint8_t>(cli.get_int("-mode", 2));
    config.strand = static_cast<int8_t>(cli.get_int("-strand", 2));
    if (config.strand != -1 && config.strand != 1 && config.strand != 2) {
        std::fprintf(stderr, "Error: -strand must be -1, 1, or 2\n");
        return 1;
    }

    // sort_score is auto-determined by mode (not a CLI option)
    switch (config.mode) {
        case 1: config.sort_score = 1; break; // stage1_score
        case 2: config.sort_score = 2; break; // chainscore
        case 3: config.sort_score = 3; break; // alnscore
        default:
            std::fprintf(stderr, "Error: -mode must be 1, 2, or 3\n");
            return 1;
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
    {
        std::string err;
        if (!parse_score_matrix(cli, stage3_config.score_matrix, err)) {
            std::fprintf(stderr, "%s\n", err.c_str());
            return 1;
        }
    }
    if (cli.has("-stage3_nthread_fetch")) {
        stage3_config.nthread_fetch = cli.get_int("-stage3_nthread_fetch", 8);
        if (stage3_config.nthread_fetch > nthread) {
            std::fprintf(stderr,
                "Error: -stage3_nthread_fetch (%d) exceeds -nthread (%d)\n",
                stage3_config.nthread_fetch, nthread);
            return 1;
        }
    } else {
        stage3_config.nthread_fetch = std::min(8, nthread);
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

    // Mode 1 consistency checks
    if (config.mode == 1) {
        // Consistency check: fractional + explicit -stage2_min_score in mode 1
        // (min_score=0 is allowed since it means adaptive)
        if (config.min_stage1_score_frac > 0 && cli.has("-stage2_min_score") &&
            config.stage2.min_score > 0) {
            std::fprintf(stderr, "Error: -stage2_min_score and fractional -stage1_min_score cannot both be specified in -mode 1\n");
            return 1;
        }

        // Consistency check: -stage2_min_score and -stage1_min_score in mode 1
        bool has_min_score = cli.has("-stage2_min_score");
        bool has_min_stage1_score = cli.has("-stage1_min_score");
        if (config.min_stage1_score_frac == 0) {
            if (has_min_score && has_min_stage1_score &&
                config.stage2.min_score != config.stage1.min_stage1_score) {
                std::fprintf(stderr, "Error: -stage2_min_score and -stage1_min_score must be the same in -mode 1\n");
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

    config.accept_qdegen = static_cast<uint8_t>(cli.get_int("-accept_qdegen", 1));

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

    if (spaced_t > 0) {
        if (!validate_spaced_seed(filter_k > 0 ? filter_k : 0, spaced_t)) {
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

    // Discover volumes. For "both" template_type we hold two contexts
    // (coding + optimal); otherwise a single context.
    const bool is_both_mode = (spaced_t > 0 && spaced_type == TemplateType::kBoth);
    std::vector<TemplateContext> ctxs;
    ctxs.reserve(2);  // pin addresses; khx_ptr depends on stable storage

    if (is_both_mode) {
        ctxs.emplace_back();
        ctxs.back().type = TemplateType::kCoding;
        ctxs.back().vol_files = discover_volumes(
            ix_prefix, filter_k, spaced_t,
            static_cast<uint8_t>(TemplateType::kCoding));
        ctxs.emplace_back();
        ctxs.back().type = TemplateType::kOptimal;
        ctxs.back().vol_files = discover_volumes(
            ix_prefix, filter_k, spaced_t,
            static_cast<uint8_t>(TemplateType::kOptimal));
        if (ctxs[0].vol_files.empty() || ctxs[1].vol_files.empty()) {
            std::fprintf(stderr,
                "Error: both-mode requires coding and optimal index files; "
                "found %zu coding, %zu optimal for prefix %s\n",
                ctxs[0].vol_files.size(), ctxs[1].vol_files.size(),
                ix_prefix.c_str());
            return 1;
        }
        if (ctxs[0].vol_files.size() != ctxs[1].vol_files.size()) {
            std::fprintf(stderr,
                "Error: coding and optimal volume counts differ (%zu vs %zu)\n",
                ctxs[0].vol_files.size(), ctxs[1].vol_files.size());
            return 1;
        }
    } else {
        ctxs.emplace_back();
        ctxs.back().type = spaced_type;
        ctxs.back().vol_files = discover_volumes(
            ix_prefix, filter_k, spaced_t,
            spaced_t > 0 ? static_cast<uint8_t>(spaced_type) : uint8_t(0));
    }

    // ctxs[0] is the canonical side (coding for both-mode, the only side otherwise).
    const auto& vol_files = ctxs[0].vol_files;

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

        // Primer pair parsing will be completed after k and seed_masks are determined.
        // For now, store the raw records for later processing.
        // We need k to be resolved first, so use a temporary flag.
        queries.reserve(primer_records.size()); // temporary; will be replaced
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

    // Helper lambda: open a set of VolumeData from discovered volume files.
    auto open_volumes = [&](const std::vector<DiscoveredVolume>& vfiles,
                            std::vector<VolumeData>& vdata,
                            bool need_kpx) -> bool {
        vdata.resize(vfiles.size());
        for (size_t vi = 0; vi < vfiles.size(); vi++) {
            const auto& vf = vfiles[vi];
            if (!vdata[vi].kix.open(vf.kix_path)) {
                std::fprintf(stderr, "Error: cannot open %s\n", vf.kix_path.c_str());
                return false;
            }
            if (need_kpx) {
                if (!vf.has_kpx) {
                    std::fprintf(stderr,
                        "Error: mode %d requires .kpx files, but %s was not found.\n"
                        "This index was built with -mode 1 (stage 1 only).\n",
                        config.mode, vf.kpx_path.c_str());
                    return false;
                }
                if (!vdata[vi].kpx.open(vf.kpx_path)) {
                    std::fprintf(stderr, "Error: cannot open %s\n", vf.kpx_path.c_str());
                    return false;
                }
            }
            if (!vdata[vi].ksx.open(vf.ksx_path)) {
                std::fprintf(stderr, "Error: cannot open %s\n", vf.ksx_path.c_str());
                return false;
            }
            vdata[vi].volume_index = vf.volume_index;
            if (filter_mode != OidFilterMode::kNone) {
                vdata[vi].filter.build(seqidlist, vdata[vi].ksx, filter_mode);
            }
        }
        return true;
    };

    const bool need_kpx = (config.mode != 1);

    // Pre-open volumes for every context.
    for (auto& ctx : ctxs) {
        if (!open_volumes(ctx.vol_files, ctx.volumes, need_kpx)) return 1;
    }

    // Open shared .khx for every context (non-fatal if missing).
    {
        auto parts = parse_index_prefix(ix_prefix);
        for (auto& ctx : ctxs) {
            ctx.khx.open(khx_path_for(parts.parent_dir, parts.db, k,
                                       spaced_t, static_cast<uint8_t>(ctx.type)));
        }
    }

    // Apply madvise budget: prioritize khx > kix dict > kpx dict > ksx.
    // Only ctxs[0]'s ksx is hinted because ksx contents are identical across
    // contexts when both sides come from the same BLAST DB.
    {
        uint64_t budget = memory_limit;
        auto try_willneed = [&budget](auto& reader) {
            uint64_t sz = reader.willneed_size();
            bool fits = (sz > 0 && budget >= sz);
            reader.apply_madvise(fits);
            if (fits) budget -= sz;
        };
        for (auto& ctx : ctxs) try_willneed(ctx.khx);
        for (auto& ctx : ctxs)
            for (auto& vd : ctx.volumes) try_willneed(vd.kix);
        for (auto& ctx : ctxs)
            for (auto& vd : ctx.volumes) try_willneed(vd.kpx);
        for (auto& vd : ctxs[0].volumes) try_willneed(vd.ksx);
        logger.info("madvise budget: %s used / %s total",
                    format_size(memory_limit - budget).c_str(),
                    format_size(memory_limit).c_str());
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

    // Phase 1: preprocess queries in parallel.  Slots are pre-allocated so
    // each task writes to ctx.pp{16,32}[qi] independently; query_pp_idx
    // collapses to the identity (qi -> qi).  Warnings and the skip-reason
    // bookkeeping run in a sequential post-pass for deterministic ordering.
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

    for (size_t qi = 0; qi < queries.size(); ++qi) {
        if (query_skip_reason[qi] != 0) {
            has_skipped = true;
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

    // Thread-local Stage1Buffer to avoid per-job allocation.
    // For "both" mode we need two buffers per thread (coding + optimal).
    uint32_t max_num_seqs = 0;
    for (const auto& vd : ctxs[0].volumes)
        max_num_seqs = std::max(max_num_seqs, vd.kix.num_sequences());

    // Determine optimal tier from actual preprocessed k-mer counts.
    uint32_t max_kmer_positions = 0;
    for (const auto& ctx : ctxs) {
        max_kmer_positions = use_uint16
            ? accumulate_max_kmer_positions(max_kmer_positions, ctx.pp16)
            : accumulate_max_kmer_positions(max_kmer_positions, ctx.pp32);
    }
    Stage1Tier tier = select_tier(max_kmer_positions, max_kmer_positions);

    auto make_tls_buf = [max_num_seqs, tier]() {
        Stage1Buffer buf;
        buf.tier = tier;
        buf.ensure_capacity(max_num_seqs);
        return buf;
    };
    tbb::enumerable_thread_specific<Stage1Buffer> tls_bufs(make_tls_buf);

    size_t num_volumes = ctxs[0].volumes.size();

    // Build VolumeBundle / QueryBundle vectors for the orchestrator. Pointers
    // alias the existing TemplateContext storage; ctxs and pp* must outlive
    // the call below.
    auto run_orchestrated = [&](auto kmer_int_tag) {
        using KmerInt = decltype(kmer_int_tag);

        std::vector<VolumeBundle<KmerInt>> volume_bundles(num_volumes);
        for (size_t vi = 0; vi < num_volumes; ++vi) {
            volume_bundles[vi].kix    = &ctxs[0].volumes[vi].kix;
            volume_bundles[vi].kpx    = &ctxs[0].volumes[vi].kpx;
            volume_bundles[vi].ksx    = &ctxs[0].volumes[vi].ksx;
            volume_bundles[vi].filter = &ctxs[0].volumes[vi].filter;
            volume_bundles[vi].volume_index = ctxs[0].volumes[vi].volume_index;
            if (is_both_mode) {
                volume_bundles[vi].kix_opt = &ctxs[1].volumes[vi].kix;
                volume_bundles[vi].kpx_opt = &ctxs[1].volumes[vi].kpx;
            }
        }

        std::vector<QueryBundle<KmerInt>> query_bundles(queries.size());
        std::vector<std::pair<size_t, size_t>> jobs;
        jobs.reserve(queries.size() * num_volumes);
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
            for (size_t vi = 0; vi < num_volumes; ++vi) {
                jobs.emplace_back(qi, vi);
            }
        }

        logger.info("Launching %zu search job(s) (Phase A) and Phase B chain extraction...",
                    jobs.size());

        tbb::task_arena arena(nthread);
        return run_search_jobs<KmerInt>(query_bundles, volume_bundles, jobs,
                                        k, config, is_both_mode, arena, tls_bufs);
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
        // Sort by query_idx (always; JSON groups by query, TSV preserves
        // input order).  When -nresult > 0 also sort within each query
        // by stage1_score descending.
        const bool need_score_sort = (config.nresult > 0);
        tbb::parallel_sort(orch_hits.begin(), orch_hits.end(),
            [need_score_sort](const OrchestratorHit& a,
                              const OrchestratorHit& b) {
                if (a.query_idx != b.query_idx)
                    return a.query_idx < b.query_idx;
                if (need_score_sort)
                    return a.cr.stage1_score > b.cr.stage1_score;
                return false;
            });

        if (need_score_sort) {
            std::vector<OrchestratorHit> truncated;
            truncated.reserve(std::min<size_t>(
                static_cast<size_t>(config.nresult) * queries.size(),
                orch_hits.size()));
            size_t cur_q = SIZE_MAX;
            uint32_t cnt = 0;
            for (auto& h : orch_hits) {
                if (h.query_idx != cur_q) { cur_q = h.query_idx; cnt = 0; }
                if (cnt < config.nresult) {
                    truncated.push_back(std::move(h));
                    ++cnt;
                }
            }
            orch_hits = std::move(truncated);
        }

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
        in.stage1_score_type = config.stage1.stage1_score_type;
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
    std::vector<OutputHit> all_hits;
    all_hits.reserve(orch_hits.size());
    for (const auto& oh_in : orch_hits) {
        const auto& cr = oh_in.cr;
        const auto& ksx_primary = ctxs[0].volumes[oh_in.volume_idx].ksx;
        OutputHit oh;
        oh.qseqid = queries[oh_in.query_idx].id;
        oh.sseqid = std::string(ksx_primary.accession(cr.seq_id));
        oh.sstrand = cr.is_reverse ? '-' : '+';
        oh.qstart = cr.q_start;
        oh.qend = cr.q_end;
        oh.sstart = cr.s_start;
        oh.send = cr.s_end;
        oh.chainscore = cr.chainscore;
        if (config.stage1.stage1_score_type == 2)
            oh.matchscore = cr.stage1_score;
        else
            oh.coverscore = cr.stage1_score;
        oh.volume = oh_in.volume_index;
        oh.oid = cr.seq_id;
        oh.qlen = static_cast<uint32_t>(queries[oh_in.query_idx].sequence.size());
        oh.slen = ksx_primary.seq_length(cr.seq_id);
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
        all_hits = run_stage3(all_hits, queries, db_path, stage3_config,
                              ctx_param.is_ratio, ctx_param.ratio, ctx_param.abs,
                              logger);
        logger.info("Stage 3 complete: %zu hits after filtering.", all_hits.size());
    }

    // Re-attach skip markers after Stage 3.
    for (auto& m : skip_markers) all_hits.push_back(std::move(m));

    // Sort and truncate final results across volumes (per query)
    if (config.nresult > 0) {
        // Sort by (query_id, sort_score desc)
        if (config.sort_score == 1) {
            tbb::parallel_sort(all_hits.begin(), all_hits.end(),
                      [](const OutputHit& a, const OutputHit& b) {
                          if (a.qseqid != b.qseqid) return a.qseqid < b.qseqid;
                          return (a.coverscore + a.matchscore) > (b.coverscore + b.matchscore);
                      });
        } else if (config.sort_score == 3) {
            tbb::parallel_sort(all_hits.begin(), all_hits.end(),
                      [](const OutputHit& a, const OutputHit& b) {
                          if (a.qseqid != b.qseqid) return a.qseqid < b.qseqid;
                          return a.alnscore > b.alnscore;
                      });
        } else {
            tbb::parallel_sort(all_hits.begin(), all_hits.end(),
                      [](const OutputHit& a, const OutputHit& b) {
                          if (a.qseqid != b.qseqid) return a.qseqid < b.qseqid;
                          return a.chainscore > b.chainscore;
                      });
        }

        // Truncate per query
        std::vector<OutputHit> truncated;
        std::string cur_qid;
        uint32_t cur_count = 0;
        for (const auto& h : all_hits) {
            if (h.qseqid != cur_qid) {
                cur_qid = h.qseqid;
                cur_count = 0;
            }
            if (cur_count < config.nresult) {
                truncated.push_back(h);
                cur_count++;
            }
        }
        all_hits = std::move(truncated);
    }
    // nresult == 0: unlimited, skip sort and truncation

    // Write output
    if (!write_all_results(output_path, all_hits, outfmt,
                           config.mode, config.stage1.stage1_score_type,
                           stage3_config.traceback, compression_level)) {
        return 1;
    }

    logger.info("Done. %zu hit(s) reported.", all_hits.size());
    return has_skipped ? 2 : 0;
}
