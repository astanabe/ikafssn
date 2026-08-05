#include "test_util.hpp"
#include "ssu_test_fixture.hpp"
#include "io/blastdb_reader.hpp"
#include "io/result_reader.hpp"
#include "io/result_writer.hpp"
#include "ikafssnretrieve/local_retriever.hpp"
#include "index/index_builder.hpp"
#include "index/kix_reader.hpp"
#include "index/kpx_reader.hpp"
#include "index/ksx_reader.hpp"
#include "search/oid_filter.hpp"
#include "search/query_preprocessor.hpp"
#include "volume_search_helper.hpp"
#include "core/config.hpp"
#include "util/logger.hpp"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

using namespace ikafssn;
using namespace ssu_fixture;

static std::string g_testdb_path;
static std::string g_test_dir;

// Runtime-extracted data
static std::string g_query_seq;
static uint32_t g_fj_oid = UINT32_MAX;

// Build index, search, write results, read results back, retrieve subsequences.
static void test_full_pipeline() {
    Stage1Buffer buf;
    std::fprintf(stderr, "-- test_full_pipeline\n");

    // Step 1: Build index
    BlastDbReader db;
    CHECK(db.open(g_testdb_path));

    Logger logger(Logger::kError);
    IndexBuilderConfig bconfig;
    bconfig.k = 7;

    std::string prefix = g_test_dir + "/ret.00.07mer";
    CHECK(build_index<uint16_t>(db, bconfig, prefix, 0, 1, "test", logger));
    db.close();

    // Step 2: Open index and search
    KixReader kix;
    KpxReader kpx;
    KsxReader ksx;
    CHECK(kix.open(prefix + ".kix"));
    CHECK(kpx.open(prefix + ".kpx"));
    CHECK(ksx.open(prefix + ".ksx"));

    OidFilter filter;
    SearchConfig config;
    config.stage1.min_stage1_score = 1;
    config.stage2.max_gap = 100;
    config.stage2.min_nhit_diag = 1;
    config.stage2.min_score = 2;

    // Query: 100bp from FJ876973.1 (extracted at runtime)
    auto qdata = preprocess_query<uint16_t>(g_query_seq, 7, nullptr, config);
    auto result = search_volume<uint16_t>(
        "query1", qdata, 7, kix, kpx, ksx, filter, config, buf);
    CHECK(!result.hits.empty());

    // Convert to OutputHit
    std::vector<OutputHit> hits;
    for (const auto& cr : result.hits) {
        OutputHit oh;
        oh.qseqid = result.query_id;
        oh.sseqid = std::string(ksx.accession(cr.seq_id));
        oh.sstrand = cr.is_reverse ? '-' : '+';
        oh.qstart = cr.q_start;
        oh.qend = cr.q_end;
        oh.sstart = cr.s_start;
        oh.send = cr.s_end;
        oh.chainscore = cr.chainscore;
        oh.volume = 0;
        hits.push_back(oh);
    }

    kix.close();
    kpx.close();
    ksx.close();

    // Step 3: Write results as tab, then read back
    std::ostringstream oss;
    write_results_tsv(oss, hits);
    std::istringstream iss(oss.str());
    bool has_volume_column = false;
    auto parsed_hits = read_results_tsv(iss, &has_volume_column);
    CHECK_EQ(parsed_hits.size(), hits.size());
    // Local retrieval routes each hit by its volume, so the column the
    // writers emit must survive the round trip.
    CHECK(has_volume_column);

    // Step 4: Retrieve subsequences from local BLAST DB
    RetrieveOptions ropts;
    ropts.context = 0;
    std::ostringstream fasta_out;
    uint32_t retrieved = retrieve_local(parsed_hits, g_testdb_path, ropts, fasta_out);

    CHECK(retrieved > 0);
    std::string fasta_str = fasta_out.str();
    CHECK(!fasta_str.empty());

    // Verify the output is valid FASTA
    CHECK(fasta_str[0] == '>');
    CHECK(fasta_str.find("query=query1") != std::string::npos);
    // The joined accession form is carried in `sseqid=`, never in the ID.
    CHECK(fasta_str.find(" sseqid=") != std::string::npos);

    // Verify extracted sequence contains only valid IUPAC nucleotide codes.
    // retrieve_local returns the DB bytes verbatim, so subjects carrying
    // degenerate bases legitimately produce IUPAC ambiguity codes.
    std::istringstream fasta_iss(fasta_str);
    std::string line;
    const std::string kIupac = "ACGTURYSWKMBDHVN";
    while (std::getline(fasta_iss, line)) {
        if (line.empty() || line[0] == '>') continue;
        for (char c : line) {
            CHECK(kIupac.find(c) != std::string::npos);
        }
    }
}

static void test_context_extension() {
    std::fprintf(stderr, "-- test_context_extension\n");

    // Create hits manually targeting ACC_FJ, hit at [100, 200)
    std::vector<OutputHit> hits;
    hits.push_back({"query1", ACC_FJ, '+', 0, 100, 100, 200, 5, 0});

    // Retrieve without context
    {
        RetrieveOptions opts;
        opts.context = 0;
        std::ostringstream out;
        uint32_t n = retrieve_local(hits, g_testdb_path, opts, out);
        CHECK_EQ(n, 1u);
        std::string fasta = out.str();
        std::istringstream iss(fasta);
        std::string hdr, seq;
        std::getline(iss, hdr);
        seq.clear();
        std::string sline;
        while (std::getline(iss, sline)) {
            if (!sline.empty() && sline[0] != '>') seq += sline;
        }
        CHECK_EQ(seq.size(), 100u); // send - sstart = 200 - 100 = 100
    }

    // Retrieve with context=10
    {
        RetrieveOptions opts;
        opts.context = 10;
        std::ostringstream out;
        uint32_t n = retrieve_local(hits, g_testdb_path, opts, out);
        CHECK_EQ(n, 1u);
        std::string fasta = out.str();
        std::istringstream iss(fasta);
        std::string hdr, seq;
        std::getline(iss, hdr);
        seq.clear();
        std::string sline;
        while (std::getline(iss, sline)) {
            if (!sline.empty() && sline[0] != '>') seq += sline;
        }
        // Range: [100-10, 200+10) = [90, 210), size = 120
        CHECK_EQ(seq.size(), 120u);
        // defline now embeds the parent accession + the
        // 1-based inclusive extracted range as `acc:start-end` (kafsss
        // style).  ext_start=90, ext_end=210 (0-based half-open) -> 91-210.
        CHECK(hdr.find(":91-210") != std::string::npos);
    }
}

static void test_context_ratio() {
    std::fprintf(stderr, "-- test_context_ratio\n");

    // Ratio mode derives the per-hit context from the hit's query length:
    // ctx = round(qlen * ratio).  qlen=100, ratio=0.1 -> ctx=10, so the range
    // [100, 200) expands to [90, 210) (size 120) — identical to abs context=10,
    // which confirms the context is computed from qlen rather than ignored.
    OutputHit h{"query1", ACC_FJ, '+', 0, 100, 100, 200, 5, 0};
    h.qlen = 100;
    std::vector<OutputHit> hits{h};

    RetrieveOptions opts;
    opts.is_ratio = true;
    opts.ratio = 0.1;
    std::ostringstream out;
    uint32_t n = retrieve_local(hits, g_testdb_path, opts, out);
    CHECK_EQ(n, 1u);

    std::string fasta = out.str();
    std::istringstream iss(fasta);
    std::string hdr, seq, sline;
    std::getline(iss, hdr);
    while (std::getline(iss, sline)) {
        if (!sline.empty() && sline[0] != '>') seq += sline;
    }
    CHECK_EQ(seq.size(), 120u);
    CHECK(hdr.find(":91-210") != std::string::npos);
}

static void test_context_clamp_start() {
    std::fprintf(stderr, "-- test_context_clamp_start\n");

    // Hit at the very start of ACC_FJ
    std::vector<OutputHit> hits;
    hits.push_back({"query1", ACC_FJ, '+', 0, 3, 0, 3, 3, 0});

    RetrieveOptions opts;
    opts.context = 10;  // context extends before position 0, should clamp
    std::ostringstream out;
    uint32_t n = retrieve_local(hits, g_testdb_path, opts, out);
    CHECK_EQ(n, 1u);
    std::string fasta = out.str();
    // defline = `>acc:1-N` after 0-based ext_start=0 is
    // converted to 1-based.  We grep for the leading ':1-' rather than the
    // upper bound (which depends on the parent length and the context
    // clamp at the upper end).
    CHECK(fasta.find(":1-") != std::string::npos);
}

static void test_reverse_strand() {
    std::fprintf(stderr, "-- test_reverse_strand\n");

    // Read first 4bp of ACC_GQ at runtime and compute expected revcomp dynamically
    BlastDbReader db;
    CHECK(db.open(g_testdb_path));
    uint32_t oid = find_oid_by_accession(db, ACC_GQ);
    CHECK(oid != UINT32_MAX);
    std::string full_seq = db.get_subsequence(oid, 0, db.seq_length(oid) - 1);
    CHECK(full_seq.size() >= 4);
    std::string first4 = full_seq.substr(0, 4);

    // Compute expected reverse complement
    std::string expected_rc;
    for (auto it = first4.rbegin(); it != first4.rend(); ++it) {
        switch (*it) {
            case 'A': expected_rc += 'T'; break;
            case 'T': expected_rc += 'A'; break;
            case 'C': expected_rc += 'G'; break;
            case 'G': expected_rc += 'C'; break;
            default:  expected_rc += 'N'; break;
        }
    }

    std::vector<OutputHit> hits;
    hits.push_back({"query1", ACC_GQ, '-', 0, 4, 0, 4, 3, 0});

    RetrieveOptions opts;
    opts.context = 0;
    std::ostringstream out;
    uint32_t n = retrieve_local(hits, g_testdb_path, opts, out);
    CHECK_EQ(n, 1u);

    std::string fasta = out.str();
    std::istringstream iss(fasta);
    std::string hdr, seq;
    std::getline(iss, hdr);
    seq.clear();
    std::string sline;
    while (std::getline(iss, sline)) {
        if (!sline.empty() && sline[0] != '>') seq += sline;
    }
    CHECK(seq == expected_rc);
    CHECK(hdr.find("strand=-") != std::string::npos);
}

static void test_missing_accession() {
    std::fprintf(stderr, "-- test_missing_accession\n");

    std::vector<OutputHit> hits;
    hits.push_back({"query1", "NONEXISTENT", '+', 0, 10, 0, 10, 5, 0});
    hits.push_back({"query1", ACC_FJ, '+', 0, 3, 0, 3, 3, 0});

    RetrieveOptions opts;
    std::ostringstream out;
    uint32_t n = retrieve_local(hits, g_testdb_path, opts, out);
    // Only seq1 should be retrieved, NONEXISTENT skipped with warning
    CHECK_EQ(n, 1u);
}

int main() {
    check_ssu_available();

    g_testdb_path = ssu_db_prefix();
    g_test_dir = "/tmp/ikafssn_retrieve_test";
    std::filesystem::create_directories(g_test_dir);

    // Extract 100bp query from FJ876973.1 at runtime
    {
        BlastDbReader db;
        CHECK(db.open(g_testdb_path));
        g_fj_oid = find_oid_by_accession(db, ACC_FJ);
        CHECK(g_fj_oid != UINT32_MAX);
        std::string full_seq = db.get_subsequence(g_fj_oid, 0, db.seq_length(g_fj_oid) - 1);
        CHECK(full_seq.size() >= 200);
        g_query_seq = full_seq.substr(100, 100);
    }

    test_full_pipeline();
    test_context_extension();
    test_context_ratio();
    test_context_clamp_start();
    test_reverse_strand();
    test_missing_accession();

    std::filesystem::remove_all(g_test_dir);

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
