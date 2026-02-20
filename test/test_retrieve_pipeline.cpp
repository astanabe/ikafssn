#include "test_util.hpp"
#include "io/blastdb_reader.hpp"
#include "io/result_reader.hpp"
#include "io/result_writer.hpp"
#include "ikafssnretrieve/local_retriever.hpp"
#include "index/index_builder.hpp"
#include "index/kix_reader.hpp"
#include "index/kpx_reader.hpp"
#include "index/ksx_reader.hpp"
#include "search/oid_filter.hpp"
#include "search/volume_searcher.hpp"
#include "core/config.hpp"
#include "util/logger.hpp"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

using namespace ikafssn;

static std::string g_testdb_path;
static std::string g_test_dir;

// Build index, search, write results, read results back, retrieve subsequences.
static void test_full_pipeline() {
    std::fprintf(stderr, "-- test_full_pipeline\n");

    // Step 1: Build index
    BlastDbReader db;
    CHECK(db.open(g_testdb_path));

    Logger logger(Logger::kError);
    IndexBuilderConfig bconfig;
    bconfig.k = 7;
    bconfig.partitions = 1;
    bconfig.buffer_size = uint64_t(1) << 30;

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
    config.stage1.max_freq = 100000;
    config.stage1.stage1_topn = 100;
    config.stage1.min_stage1_score = 1;
    config.stage2.max_gap = 100;
    config.stage2.min_diag_hits = 1;
    config.stage2.min_score = 2;
    config.num_results = 50;

    // Query: first 24 bases of seq1
    std::string query_seq = "ACGTACGTACGTACGTACGTACGT";
    auto result = search_volume<uint16_t>(
        "query1", query_seq, 7, kix, kpx, ksx, filter, config);
    CHECK(!result.hits.empty());

    // Convert to OutputHit
    std::vector<OutputHit> hits;
    for (const auto& cr : result.hits) {
        OutputHit oh;
        oh.query_id = result.query_id;
        oh.accession = std::string(ksx.accession(cr.seq_id));
        oh.strand = cr.is_reverse ? '-' : '+';
        oh.q_start = cr.q_start;
        oh.q_end = cr.q_end;
        oh.s_start = cr.s_start;
        oh.s_end = cr.s_end;
        oh.score = cr.score;
        oh.volume = 0;
        hits.push_back(oh);
    }

    kix.close();
    kpx.close();
    ksx.close();

    // Step 3: Write results as tab, then read back
    std::ostringstream oss;
    write_results_tab(oss, hits);
    std::istringstream iss(oss.str());
    auto parsed_hits = read_results_tab(iss);
    CHECK_EQ(parsed_hits.size(), hits.size());

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

    // Verify extracted sequence contains only valid bases
    std::istringstream fasta_iss(fasta_str);
    std::string line;
    while (std::getline(fasta_iss, line)) {
        if (line.empty() || line[0] == '>') continue;
        for (char c : line) {
            CHECK(c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N');
        }
    }
}

static void test_context_extension() {
    std::fprintf(stderr, "-- test_context_extension\n");

    // Create hits manually targeting seq1 (ACGTACGTACGTACGTACGTACGTACGTACGT, 32bp)
    // Match region: [4, 11] (8 bases)
    std::vector<OutputHit> hits;
    hits.push_back({"query1", "seq1", '+', 0, 7, 4, 11, 5, 0});

    // Retrieve without context
    {
        RetrieveOptions opts;
        opts.context = 0;
        std::ostringstream out;
        uint32_t n = retrieve_local(hits, g_testdb_path, opts, out);
        CHECK_EQ(n, 1u);
        std::string fasta = out.str();
        // Extract sequence (skip header line)
        std::istringstream iss(fasta);
        std::string hdr, seq;
        std::getline(iss, hdr);
        seq.clear();
        std::string sline;
        while (std::getline(iss, sline)) {
            if (!sline.empty() && sline[0] != '>') seq += sline;
        }
        CHECK_EQ(seq.size(), 8u); // s_end - s_start + 1 = 11 - 4 + 1 = 8
    }

    // Retrieve with context=3
    {
        RetrieveOptions opts;
        opts.context = 3;
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
        // Range: [4-3, 11+3] = [1, 14], size = 14
        CHECK_EQ(seq.size(), 14u);
        CHECK(hdr.find("range=1-14") != std::string::npos);
    }
}

static void test_context_clamp_start() {
    std::fprintf(stderr, "-- test_context_clamp_start\n");

    // Hit at the very start of the sequence
    std::vector<OutputHit> hits;
    hits.push_back({"query1", "seq1", '+', 0, 3, 0, 3, 3, 0});

    RetrieveOptions opts;
    opts.context = 10;  // context extends before position 0, should clamp
    std::ostringstream out;
    uint32_t n = retrieve_local(hits, g_testdb_path, opts, out);
    CHECK_EQ(n, 1u);
    std::string fasta = out.str();
    CHECK(fasta.find("range=0-") != std::string::npos);
}

static void test_reverse_strand() {
    std::fprintf(stderr, "-- test_reverse_strand\n");

    // seq1 = "ACGTACGTACGTACGTACGTACGTACGTACGT" (32bp)
    // Take a region [0,3] -> "ACGT", reverse complement = "ACGT" (palindrome)
    // Take [0,7] -> "ACGTACGT", revcomp = "ACGTACGT" (also palindrome)
    // Use a non-palindromic region from seq5:
    // seq5 = "ACGTACGTAACCGGTTAACCGGTTACGTACGTAACCGGTTAACCGGTT" (48bp)
    // Region [8,15] = "AACCGGTT", revcomp = "AACCGGTT" (palindromic too!)
    // Region [0,3] of seq4 = "ACAC", revcomp = "GTGT"
    std::vector<OutputHit> hits;
    hits.push_back({"query1", "seq4", '-', 0, 3, 0, 3, 3, 0});

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
    // seq4[0..3] = "ACAC", revcomp = "GTGT"
    CHECK(seq == "GTGT");
    CHECK(hdr.find("strand=-") != std::string::npos);
}

static void test_missing_accession() {
    std::fprintf(stderr, "-- test_missing_accession\n");

    std::vector<OutputHit> hits;
    hits.push_back({"query1", "NONEXISTENT", '+', 0, 10, 0, 10, 5, 0});
    hits.push_back({"query1", "seq1", '+', 0, 3, 0, 3, 3, 0});

    RetrieveOptions opts;
    std::ostringstream out;
    uint32_t n = retrieve_local(hits, g_testdb_path, opts, out);
    // Only seq1 should be retrieved, NONEXISTENT skipped with warning
    CHECK_EQ(n, 1u);
}

int main() {
    g_testdb_path = std::string(SOURCE_DIR) + "/test/testdata/testdb";
    g_test_dir = "/tmp/ikafssn_retrieve_test";
    std::filesystem::create_directories(g_test_dir);

    test_full_pipeline();
    test_context_extension();
    test_context_clamp_start();
    test_reverse_strand();
    test_missing_accession();

    std::filesystem::remove_all(g_test_dir);

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
