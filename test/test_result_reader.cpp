#include "test_util.hpp"
#include "io/output_coords.hpp"
#include "io/result_reader.hpp"
#include "io/result_writer.hpp"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

using namespace ikafssn;

static std::string g_test_dir;

static void write_file(const std::string& path, const std::string& content) {
    std::ofstream f(path);
    f << content;
}

// Rows carry the BLAST `-outfmt 6` convention (1-based inclusive,
// query-relative, sstart > send on the reverse strand); the reader hands
// back the internal 0-based half-open, query-strand-relative form.
static void test_basic_parse() {
    std::fprintf(stderr, "-- test_basic_parse\n");

    std::string path = g_test_dir + "/basic.tsv";
    write_file(path,
        "# qseqid\tsseqid\tsstrand\tqstart\tqend\tqlen\tsstart\tsend\tslen\tcoverscore\tchainscore\tvolume\n"
        "query1\tACC001\t+\t1\t49\t500\t101\t149\t2000\t5\t15\t0\n"
        "query1\tACC002\t-\t11\t40\t500\t229\t201\t3000\t3\t10\t1\n"
    );

    auto results = read_results_tsv(path);
    CHECK_EQ(results.size(), 2u);

    CHECK(results[0].qseqid == "query1");
    CHECK(results[0].sseqid == "ACC001");
    CHECK_EQ(results[0].sstrand, '+');
    CHECK_EQ(results[0].qstart, 0u);
    CHECK_EQ(results[0].qend, 49u);
    CHECK_EQ(results[0].qlen, 500u);
    CHECK_EQ(results[0].sstart, 100u);
    CHECK_EQ(results[0].send, 149u);
    CHECK_EQ(results[0].slen, 2000u);
    CHECK_EQ(results[0].coverscore, 5u);
    CHECK_EQ(results[0].chainscore, 15u);
    CHECK_EQ(results[0].volume, 0u);

    CHECK(results[1].qseqid == "query1");
    CHECK(results[1].sseqid == "ACC002");
    CHECK_EQ(results[1].sstrand, '-');
    // Query interval [11, 40] of a 500 bp query folds to [460, 490) on the
    // reverse complement; the subject interval un-swaps to [200, 229).
    CHECK_EQ(results[1].qstart, 460u);
    CHECK_EQ(results[1].qend, 490u);
    CHECK_EQ(results[1].qlen, 500u);
    CHECK_EQ(results[1].sstart, 200u);
    CHECK_EQ(results[1].send, 229u);
    CHECK_EQ(results[1].slen, 3000u);
    CHECK_EQ(results[1].coverscore, 3u);
    CHECK_EQ(results[1].chainscore, 10u);
    CHECK_EQ(results[1].volume, 1u);
}

static void test_skip_header_and_blank() {
    std::fprintf(stderr, "-- test_skip_header_and_blank\n");

    std::string path = g_test_dir + "/header.tsv";
    write_file(path,
        "# comment line\n"
        "# qseqid\tsseqid\tsstrand\tqstart\tqend\tqlen\tsstart\tsend\tslen\tcoverscore\tchainscore\tvolume\n"
        "\n"
        "query1\tACC001\t+\t1\t49\t500\t101\t149\t2000\t5\t15\t0\n"
        "\n"
    );

    auto results = read_results_tsv(path);
    CHECK_EQ(results.size(), 1u);
    CHECK(results[0].sseqid == "ACC001");
}

static void test_invalid_lines() {
    std::fprintf(stderr, "-- test_invalid_lines\n");

    std::string path = g_test_dir + "/invalid.tsv";
    write_file(path,
        "# qseqid\tsseqid\tsstrand\tqstart\tqend\tqlen\tsstart\tsend\tslen\tcoverscore\tchainscore\tvolume\n"
        "query1\tACC001\t+\t1\t49\t500\t101\t149\t2000\t5\t15\t0\n"
        "too_few_fields\tACC002\n"
        "query2\tACC003\tX\t1\t49\t500\t101\t149\t2000\t5\t15\t0\n"   // bad strand
        "query3\tACC004\t+\tabc\t49\t500\t101\t149\t2000\t5\t15\t0\n"  // bad number
        "query4\tACC005\t-\t6\t55\t500\t350\t301\t2000\t8\t20\t2\n"
    );

    auto results = read_results_tsv(path);
    CHECK_EQ(results.size(), 2u);
    CHECK(results[0].sseqid == "ACC001");
    CHECK(results[1].sseqid == "ACC005");
}

static void test_empty_input() {
    std::fprintf(stderr, "-- test_empty_input\n");

    std::string path = g_test_dir + "/empty.tsv";
    write_file(path, "");

    auto results = read_results_tsv(path);
    CHECK_EQ(results.size(), 0u);
}

static void test_header_only() {
    std::fprintf(stderr, "-- test_header_only\n");

    std::string path = g_test_dir + "/header_only.tsv";
    write_file(path,
        "# qseqid\tsseqid\tsstrand\tqstart\tqend\tsstart\tsend\tscore\tvolume\n"
    );

    auto results = read_results_tsv(path);
    CHECK_EQ(results.size(), 0u);
}

static void test_stream_interface() {
    std::fprintf(stderr, "-- test_stream_interface\n");

    std::istringstream iss(
        "# qseqid\tsseqid\tsstrand\tqstart\tqend\tqlen\tsstart\tsend\tslen\tcoverscore\tchainscore\tvolume\n"
        "q1\tA1\t+\t1\t10\t100\t21\t30\t200\t3\t5\t0\n"
        "q2\tA2\t-\t6\t15\t100\t35\t26\t200\t4\t8\t1\n"
    );

    auto results = read_results_tsv(iss);
    CHECK_EQ(results.size(), 2u);
    CHECK(results[0].qseqid == "q1");
    CHECK(results[1].qseqid == "q2");
}

static void test_roundtrip() {
    std::fprintf(stderr, "-- test_roundtrip\n");

    // Create hits, write them, then read back
    std::vector<OutputHit> hits;
    hits.push_back({"queryA", "ACC100", '+', 0, 99, 500, 599, 25, 8, 0});
    hits.push_back({"queryA", "ACC200", '-', 10, 89, 1000, 1079, 18, 5, 1});
    hits.push_back({"queryB", "ACC300", '+', 0, 49, 0, 49, 12, 3, 0});

    std::ostringstream oss;
    write_results_tsv(oss, hits);

    std::istringstream iss(oss.str());
    auto read_back = read_results_tsv(iss);

    CHECK_EQ(read_back.size(), 3u);
    for (size_t i = 0; i < 3; i++) {
        CHECK(read_back[i].qseqid == hits[i].qseqid);
        CHECK(read_back[i].sseqid == hits[i].sseqid);
        CHECK_EQ(read_back[i].sstrand, hits[i].sstrand);
        CHECK_EQ(read_back[i].qstart, hits[i].qstart);
        CHECK_EQ(read_back[i].qend, hits[i].qend);
        CHECK_EQ(read_back[i].sstart, hits[i].sstart);
        CHECK_EQ(read_back[i].send, hits[i].send);
        CHECK_EQ(read_back[i].coverscore, hits[i].coverscore);
        CHECK_EQ(read_back[i].chainscore, hits[i].chainscore);
        CHECK_EQ(read_back[i].volume, hits[i].volume);
    }
}

static void test_roundtrip_mode3_no_traceback() {
    std::fprintf(stderr, "-- test_roundtrip_mode3_no_traceback\n");

    // Mode 3, no traceback: qstart/sstart should NOT be in output
    std::vector<OutputHit> hits;
    OutputHit h;
    h.qseqid = "qryM3";
    h.sseqid = "ACC_M3";
    h.sstrand = '+';
    h.qstart = 5;    // will NOT be written
    h.qend = 95;
    h.qlen = 100;
    h.sstart = 200;  // will NOT be written
    h.send = 800;
    h.slen = 5000;
    h.coverscore = 10;
    h.chainscore = 50;
    h.alnscore = 120;
    h.volume = 2;
    hits.push_back(h);

    std::ostringstream oss;
    write_results_tsv(oss, hits, /*mode=*/3, /*stage3_traceback=*/false);

    // Verify qstart/sstart are not in header
    std::string output = oss.str();
    CHECK(output.find("qstart") == std::string::npos);
    CHECK(output.find("sstart") == std::string::npos);
    CHECK(output.find("qend") != std::string::npos);
    CHECK(output.find("send") != std::string::npos);

    // Read back
    std::istringstream iss(output);
    auto read_back = read_results_tsv(iss);
    CHECK_EQ(read_back.size(), 1u);
    CHECK(read_back[0].qseqid == "qryM3");
    CHECK(read_back[0].sseqid == "ACC_M3");
    CHECK_EQ(read_back[0].sstrand, '+');
    CHECK_EQ(read_back[0].qstart, 0u);   // not present -> default 0
    CHECK_EQ(read_back[0].qend, 95u);
    CHECK_EQ(read_back[0].qlen, 100u);
    CHECK_EQ(read_back[0].sstart, 0u);   // not present -> default 0
    CHECK_EQ(read_back[0].send, 800u);
    CHECK_EQ(read_back[0].slen, 5000u);
    CHECK_EQ(read_back[0].coverscore, 10u);
    CHECK_EQ(read_back[0].chainscore, 50u);
    CHECK_EQ(read_back[0].alnscore, 120);
    CHECK_EQ(read_back[0].volume, 2u);
}

static void test_roundtrip_mode3_traceback() {
    std::fprintf(stderr, "-- test_roundtrip_mode3_traceback\n");

    // Mode 3, with traceback: qstart/sstart SHOULD be present
    std::vector<OutputHit> hits;
    OutputHit h;
    h.qseqid = "qryTB";
    h.sseqid = "ACC_TB";
    h.sstrand = '-';
    h.qstart = 3;
    h.qend = 97;
    h.qlen = 100;
    h.sstart = 150;
    h.send = 750;
    h.slen = 3000;
    h.coverscore = 7;
    h.chainscore = 40;
    h.alnscore = 200;
    h.ppositive = 95.5;
    h.npositive = 90;
    h.nnegative = 4;
    h.cigar = "50M2I48M";
    h.qseq = "ACGT";
    h.sseq = "ACGT";
    h.volume = 1;
    hits.push_back(h);

    std::ostringstream oss;
    write_results_tsv(oss, hits, /*mode=*/3, /*stage3_traceback=*/true);

    std::string output = oss.str();
    CHECK(output.find("qstart") != std::string::npos);
    CHECK(output.find("sstart") != std::string::npos);

    std::istringstream iss(output);
    auto read_back = read_results_tsv(iss);
    CHECK_EQ(read_back.size(), 1u);
    CHECK_EQ(read_back[0].qstart, 3u);
    CHECK_EQ(read_back[0].qend, 97u);
    CHECK_EQ(read_back[0].sstart, 150u);
    CHECK_EQ(read_back[0].send, 750u);
    CHECK_EQ(read_back[0].alnscore, 200);
    CHECK(read_back[0].cigar == "50M2I48M");
    CHECK(read_back[0].qseq == "ACGT");
    CHECK(read_back[0].sseq == "ACGT");
}

// The exact bytes of the two fields whose formatting is not a plain
// unsigned decimal: ppositive (six significant digits, C locale) and
// alnscore (signed).
static void test_numeric_field_formatting() {
    std::fprintf(stderr, "-- test_numeric_field_formatting\n");

    struct Case { double ppositive; int32_t alnscore; const char* pp; const char* aln; };
    const Case cases[] = {
        {0.0,                 0,           "0",           "0"},
        {100.0,               1000,        "100",         "1000"},
        {94.73684210526316,  -5,           "94.7368",     "-5"},
        {33.333333333333336, -2147483648,  "33.3333",     "-2147483648"},
        {99.9999994,          2147483647,  "100",         "2147483647"},
        {66.66666666666667,  -1,           "66.6667",     "-1"},
        {0.000123456789,      7,           "0.000123457", "7"},
        {1e-5,               -12345,       "1e-05",       "-12345"},
    };

    for (const auto& c : cases) {
        std::vector<OutputHit> hits;
        OutputHit h;
        h.qseqid = "q";
        h.sseqid = "s";
        h.sstrand = '+';
        h.qstart = 0; h.qend = 2; h.qlen = 3;
        h.sstart = 3; h.send = 5; h.slen = 6;
        h.coverscore = 7; h.chainscore = 8;
        h.alnscore = c.alnscore;
        h.ppositive = c.ppositive;
        h.npositive = 9; h.nnegative = 10;
        h.cigar = "1M"; h.qseq = "A"; h.sseq = "A";
        h.volume = 0;
        hits.push_back(h);

        std::ostringstream tsv;
        write_results_tsv(tsv, hits, /*mode=*/3, /*stage3_traceback=*/true);
        std::string want_tsv = std::string("q\ts\t+\t1\t2\t3\t4\t5\t6\t7\t8\t") +
                               c.aln + "\t" + c.pp + "\t9\t10\t1M\tA\tA\t0\n";
        CHECK(tsv.str().size() > want_tsv.size());
        CHECK(tsv.str().compare(tsv.str().size() - want_tsv.size(),
                                want_tsv.size(), want_tsv) == 0);

        std::ostringstream json;
        write_results_json(json, hits, /*mode=*/3, /*stage3_traceback=*/true);
        CHECK(json.str().find(std::string("\"alnscore\": ") + c.aln + ",") !=
              std::string::npos);
        CHECK(json.str().find(std::string("\"ppositive\": ") + c.pp + ",") !=
              std::string::npos);
    }
}

// JSON string fields carry whatever the FASTA header or a failure reason
// held, so the escaping has to survive a round trip through the writer.
static void test_json_string_escaping() {
    std::fprintf(stderr, "-- test_json_string_escaping\n");

    std::vector<OutputHit> hits;
    OutputHit h;
    h.qseqid = "q\"quoted\"";
    h.sseqid = "acc\\back";
    h.sstrand = '+';
    h.qstart = 1; h.qend = 2; h.qlen = 3;
    h.sstart = 4; h.send = 5; h.slen = 6;
    h.coverscore = 7; h.chainscore = 8;
    h.alnscore = 9; h.ppositive = 50.0;
    h.npositive = 1; h.nnegative = 1;
    h.cigar = "a\tb\nc\rd";
    h.qseq = "A"; h.sseq = "A";
    h.volume = 0;
    hits.push_back(h);

    OutputHit sk;
    sk.qseqid = "skipped\tquery";
    sk.qlen = 10;
    sk.skip_reason = ikafssn::kSkipQueryTooShort;
    sk.skip_detail = "line1\nline2 \"x\"";
    hits.push_back(sk);

    std::ostringstream oss;
    write_results_json(oss, hits, /*mode=*/3, /*stage3_traceback=*/true);
    const std::string out = oss.str();

    CHECK(out.find("\"qseqid\": \"q\\\"quoted\\\"\"") != std::string::npos);
    CHECK(out.find("\"sseqid\": \"acc\\\\back\"") != std::string::npos);
    CHECK(out.find("\"cigar\": \"a\\tb\\nc\\rd\"") != std::string::npos);
    CHECK(out.find("\"qseqid\": \"skipped\\tquery\"") != std::string::npos);
    CHECK(out.find("\"skip_detail\": \"line1\\nline2 \\\"x\\\"\"") != std::string::npos);
    // The raw control characters must not survive into the JSON.
    CHECK(out.find('\t') == std::string::npos);
    CHECK(out.find('\r') == std::string::npos);
}

// `nthread` must never change the output, whichever path the writer takes.
static void test_writer_thread_count_invariance() {
    std::fprintf(stderr, "-- test_writer_thread_count_invariance\n");

    // Uneven hits per query, including a query that is only a skip
    // sentinel.  Small enough that the writers stay on the serial path.
    std::vector<OutputHit> hits;
    for (int q = 0; q < 40; q++) {
        std::string qid = "query" + std::to_string(q);
        if (q % 11 == 3) {   // skipped query: no hits, sentinel row only
            OutputHit s;
            s.qseqid = qid;
            s.qlen = 100 + q;
            s.skip_reason = ikafssn::kSkipQueryTooShort;
            s.skip_detail = "too short\tand \"quoted\"";
            hits.push_back(s);
            continue;
        }
        for (int i = 0; i < (q % 7) * 900 + 1; i++) {
            OutputHit h;
            h.qseqid = qid;
            h.sseqid = "ACC" + std::to_string(i);
            h.sstrand = (i % 2) ? '-' : '+';
            h.qstart = i; h.qend = i + 40; h.qlen = 100 + q;
            h.sstart = i * 2; h.send = i * 2 + 40; h.slen = 5000 + i;
            h.coverscore = i % 31; h.chainscore = i % 17;
            h.alnscore = (i % 5 == 0) ? -(i + 1) : i;
            h.ppositive = 100.0 * (i % 19) / 19.0;
            h.npositive = i % 19; h.nnegative = 19 - (i % 19);
            h.cigar = "40M"; h.qseq = "ACGT"; h.sseq = "ACGA";
            h.volume = static_cast<uint16_t>(i % 3);
            hits.push_back(h);
        }
    }

    struct Case { uint8_t mode; bool tb; };
    const Case cases[] = {{1, false}, {2, false}, {3, false}, {3, true}};
    for (const auto& c : cases) {
        for (int nthread : {2, 8}) {
            std::ostringstream a, b;
            write_results_tsv(a, hits, c.mode, c.tb, /*nthread=*/1);
            write_results_tsv(b, hits, c.mode, c.tb, nthread);
            CHECK(a.str() == b.str());

            std::ostringstream ja, jb;
            write_results_json(ja, hits, c.mode, c.tb, /*nthread=*/1);
            write_results_json(jb, hits, c.mode, c.tb, nthread);
            CHECK(ja.str() == jb.str());

            std::ostringstream fa, fb;
            write_results_json_fragment(fa, hits, c.mode, c.tb, /*nthread=*/1);
            write_results_json_fragment(fb, hits, c.mode, c.tb, nthread);
            CHECK(fa.str() == fb.str());
        }
    }
}

// Above kParallelMinWeight the writers chunk the work across threads and
// format it a wave at a time, so chunk and wave boundaries have to fall
// inside the output without disturbing it.
static void test_writer_parallel_chunking() {
    std::fprintf(stderr, "-- test_writer_parallel_chunking\n");

    std::vector<OutputHit> hits;
    for (int q = 0; q < 40; q++) {
        std::string qid = "query" + std::to_string(q);
        if (q % 11 == 3) {
            OutputHit s;
            s.qseqid = qid;
            s.qlen = 100 + q;
            s.skip_reason = ikafssn::kSkipQueryTooShort;
            s.skip_detail = "too short";
            hits.push_back(s);
            continue;
        }
        // Uneven per query so a chunk boundary lands inside one, and large
        // enough in total to cross both the parallel and the wave threshold.
        for (int i = 0; i < (q % 7) * 2600 + 1; i++) {
            OutputHit h;
            h.qseqid = qid;
            h.sseqid = "ACC" + std::to_string(i);
            h.sstrand = (i % 2) ? '-' : '+';
            h.qstart = i; h.qend = i + 40; h.qlen = 100 + q;
            h.sstart = i * 2; h.send = i * 2 + 40; h.slen = 5000 + i;
            h.coverscore = i % 31; h.chainscore = i % 17;
            h.volume = static_cast<uint16_t>(i % 3);
            hits.push_back(h);
        }
    }

    {
        std::ostringstream a, b;
        write_results_tsv(a, hits, /*mode=*/2, /*stage3_traceback=*/false, 1);
        write_results_tsv(b, hits, /*mode=*/2, /*stage3_traceback=*/false, 8);
        CHECK(a.str() == b.str());
    }
    {
        std::ostringstream a, b;
        write_results_json(a, hits, /*mode=*/2, /*stage3_traceback=*/false, 1);
        write_results_json(b, hits, /*mode=*/2, /*stage3_traceback=*/false, 8);
        CHECK(a.str() == b.str());
    }
}

static void test_header_reordered_columns() {
    std::fprintf(stderr, "-- test_header_reordered_columns\n");

    // Header with columns in a different order
    std::string path = g_test_dir + "/reordered.tsv";
    write_file(path,
        "# sseqid\tqseqid\tvolume\tsstrand\tcoverscore\tqlen\tslen\n"
        "ACC_R1\tqR1\t3\t+\t12\t400\t5000\n"
        "ACC_R2\tqR2\t0\t-\t8\t300\t4000\n"
    );

    auto results = read_results_tsv(path);
    CHECK_EQ(results.size(), 2u);

    CHECK(results[0].qseqid == "qR1");
    CHECK(results[0].sseqid == "ACC_R1");
    CHECK_EQ(results[0].sstrand, '+');
    CHECK_EQ(results[0].coverscore, 12u);
    CHECK_EQ(results[0].qlen, 400u);
    CHECK_EQ(results[0].slen, 5000u);
    CHECK_EQ(results[0].volume, 3u);
    // Missing columns default to 0
    CHECK_EQ(results[0].qstart, 0u);
    CHECK_EQ(results[0].sstart, 0u);
    CHECK_EQ(results[0].chainscore, 0u);

    CHECK(results[1].qseqid == "qR2");
    CHECK(results[1].sseqid == "ACC_R2");
    CHECK_EQ(results[1].sstrand, '-');
    CHECK_EQ(results[1].coverscore, 8u);
    CHECK_EQ(results[1].volume, 0u);
}


static void test_no_header_fallback() {
    std::fprintf(stderr, "-- test_no_header_fallback\n");

    // No valid header line — should fall back to the field-count parser,
    // which applies the same BLAST-convention inverse as the mapped path.
    std::istringstream iss(
        "q1\tA1\t+\t1\t10\t100\t21\t30\t200\t3\t5\t0\n"
        "q2\tA2\t-\t6\t15\t100\t35\t26\t200\t4\t8\t1\n"
    );

    auto results = read_results_tsv(iss);
    CHECK_EQ(results.size(), 2u);
    CHECK(results[0].qseqid == "q1");
    CHECK_EQ(results[0].qstart, 0u);
    CHECK_EQ(results[0].qend, 10u);
    CHECK_EQ(results[0].sstart, 20u);
    CHECK_EQ(results[0].send, 30u);
    CHECK(results[1].qseqid == "q2");
    CHECK_EQ(results[1].qstart, 85u);
    CHECK_EQ(results[1].qend, 95u);
    CHECK_EQ(results[1].sstart, 25u);
    CHECK_EQ(results[1].send, 35u);
}

// A missing volume column reads back as volume 0, so `has_volume_column` is
// the only way ikafssnretrieve can tell the two apart.
static void test_has_volume_column() {
    std::fprintf(stderr, "-- test_has_volume_column\n");

    const char* const kHeader =
        "# qseqid\tsseqid\tsstrand\tqstart\tqend\tqlen\tsstart\tsend\tslen"
        "\tcoverscore\tchainscore";

    {
        std::istringstream iss(
            std::string(kHeader) + "\tvolume\n"
            "q1\tA1\t+\t1\t10\t100\t21\t30\t200\t3\t5\t7\n");
        bool has_volume = false;
        auto results = read_results_tsv(iss, &has_volume);
        CHECK_EQ(results.size(), 1u);
        CHECK(has_volume);
        CHECK_EQ(results[0].volume, 7u);
    }

    {
        std::istringstream iss(
            std::string(kHeader) + "\n"
            "q1\tA1\t+\t1\t10\t100\t21\t30\t200\t3\t5\n");
        bool has_volume = true;
        auto results = read_results_tsv(iss, &has_volume);
        CHECK_EQ(results.size(), 1u);
        CHECK(!has_volume);
        CHECK_EQ(results[0].volume, 0u);
    }

    // The field-count fallback recognises no layout without a volume field.
    {
        std::istringstream iss("q1\tA1\t+\t1\t10\t100\t21\t30\t200\t3\t5\t7\n");
        bool has_volume = false;
        auto results = read_results_tsv(iss, &has_volume);
        CHECK_EQ(results.size(), 1u);
        CHECK(has_volume);
        CHECK_EQ(results[0].volume, 7u);
    }

    // Nothing parsed leaves it false.
    {
        std::istringstream iss("");
        bool has_volume = true;
        auto results = read_results_tsv(iss, &has_volume);
        CHECK_EQ(results.size(), 0u);
        CHECK(!has_volume);
    }
}

static void test_output_coords_roundtrip() {
    std::fprintf(stderr, "-- test_output_coords_roundtrip\n");

    struct Case { uint32_t q_start, q_end, s_start, s_end, qlen; };
    const Case cases[] = {
        {  0, 200,   0, 300, 200},   // both intervals start at the origin
        { 10, 200,  55, 300, 200},   // q_end == qlen
        {  0,  17,   0,  17,  17},   // whole query, whole subject prefix
        { 41,  73, 900, 932, 512},
        {  0,   1,   0,   1,   1},   // single-base query
    };

    for (const auto& c : cases) {
        for (bool is_reverse : {false, true}) {
            OutputCoords o = to_output_coords(c.q_start, c.q_end,
                                              c.s_start, c.s_end,
                                              c.qlen, is_reverse);
            // BLAST always reports ascending query coordinates.
            CHECK(o.qstart <= o.qend);
            CHECK(o.qstart >= 1);
            CHECK(o.sstart >= 1 && o.send >= 1);
            // The subject runs in alignment order.
            if (is_reverse) CHECK(o.sstart >= o.send);
            else            CHECK(o.sstart <= o.send);

            uint32_t q_start = 0, q_end = 0, s_start = 0, s_end = 0;
            to_internal_coords(o.qstart, o.qend, o.sstart, o.send,
                               c.qlen, is_reverse,
                               q_start, q_end, s_start, s_end);
            CHECK_EQ(q_start, c.q_start);
            CHECK_EQ(q_end, c.q_end);
            CHECK_EQ(s_start, c.s_start);
            CHECK_EQ(s_end, c.s_end);
        }
    }
}

static void test_result_reader_blast_coords() {
    std::fprintf(stderr, "-- test_result_reader_blast_coords\n");

    std::vector<OutputHit> hits;
    OutputHit fwd;
    fwd.qseqid = "qF"; fwd.sseqid = "ACC_F"; fwd.sstrand = '+';
    fwd.qstart = 0; fwd.qend = 200; fwd.qlen = 200;
    fwd.sstart = 100; fwd.send = 300; fwd.slen = 1725;
    fwd.coverscore = 190; fwd.chainscore = 190; fwd.volume = 0;
    hits.push_back(fwd);

    OutputHit rev = fwd;
    rev.qseqid = "qR"; rev.sstrand = '-';
    hits.push_back(rev);

    std::ostringstream tsv;
    write_results_tsv(tsv, hits);

    // The exact BLAST-convention rows, matching what blastn -outfmt 6 emits
    // for the same alignment.
    CHECK(tsv.str().find("qF\tACC_F\t+\t1\t200\t200\t101\t300\t1725\t") !=
          std::string::npos);
    CHECK(tsv.str().find("qR\tACC_F\t-\t1\t200\t200\t300\t101\t1725\t") !=
          std::string::npos);

    // JSON carries the same coordinates.  There is no JSON reader, so the
    // written text is what gets checked.
    std::ostringstream json;
    write_results_json(json, hits);
    CHECK(json.str().find("\"sstart\": 101") != std::string::npos);
    CHECK(json.str().find("\"send\": 300") != std::string::npos);
    CHECK(json.str().find("\"sstart\": 300") != std::string::npos);
    CHECK(json.str().find("\"send\": 101") != std::string::npos);
    CHECK(json.str().find("\"qstart\": 1") != std::string::npos);
    CHECK(json.str().find("\"qend\": 200") != std::string::npos);

    std::istringstream iss(tsv.str());
    auto read_back = read_results_tsv(iss);
    CHECK_EQ(read_back.size(), 2u);
    for (size_t i = 0; i < read_back.size() && i < hits.size(); i++) {
        CHECK_EQ(read_back[i].sstrand, hits[i].sstrand);
        CHECK_EQ(read_back[i].qstart, hits[i].qstart);
        CHECK_EQ(read_back[i].qend, hits[i].qend);
        CHECK_EQ(read_back[i].sstart, hits[i].sstart);
        CHECK_EQ(read_back[i].send, hits[i].send);
    }
}

static void test_windows_line_endings() {
    std::fprintf(stderr, "-- test_windows_line_endings\n");

    std::string path = g_test_dir + "/crlf.tsv";
    write_file(path,
        "# qseqid\tsseqid\tsstrand\tqstart\tqend\tqlen\tsstart\tsend\tslen\tcoverscore\tchainscore\tvolume\r\n"
        "q1\tA1\t+\t1\t10\t100\t21\t30\t200\t3\t5\t0\r\n"
    );

    auto results = read_results_tsv(path);
    CHECK_EQ(results.size(), 1u);
    CHECK(results[0].qseqid == "q1");
    CHECK(results[0].sseqid == "A1");
}

// The JSON writer emits one object per qseqid, in the order the queries first
// appear, and each object carries that query's hits in input order — even when
// the hits of one query are not contiguous.
static void test_json_query_first_appearance_order() {
    std::fprintf(stderr, "-- test_json_query_first_appearance_order\n");

    auto make = [](const char* q, const char* s) {
        OutputHit h;
        h.qseqid = q; h.sseqid = s; h.sstrand = '+';
        h.qstart = 0; h.qend = 10; h.qlen = 10;
        h.sstart = 0; h.send = 10; h.slen = 100;
        h.coverscore = 1; h.chainscore = 2; h.volume = 0;
        return h;
    };
    std::vector<OutputHit> hits{make("qB", "sB1"), make("qA", "sA1"),
                                make("qB", "sB2")};

    std::ostringstream oss;
    write_results_json(oss, hits, /*mode=*/2, /*stage3_traceback=*/false);
    const std::string out = oss.str();
    const std::size_t npos = std::string::npos;

    const std::size_t p_qb = out.find("\"qseqid\": \"qB\"");
    const std::size_t p_qa = out.find("\"qseqid\": \"qA\"");
    CHECK(p_qb != npos);
    CHECK(p_qa != npos);
    if (p_qb == npos || p_qa == npos) return;
    // Two objects, qB's first: the lexicographic order would put qA first.
    CHECK(p_qb < p_qa);
    CHECK(out.find("\"qseqid\": \"qB\"", p_qb + 1) == npos);
    CHECK(out.find("\"qseqid\": \"qA\"", p_qa + 1) == npos);

    const std::size_t p_sb1 = out.find("\"sseqid\": \"sB1\"");
    const std::size_t p_sb2 = out.find("\"sseqid\": \"sB2\"");
    const std::size_t p_sa1 = out.find("\"sseqid\": \"sA1\"");
    CHECK(p_sb1 != npos);
    CHECK(p_sb2 != npos);
    CHECK(p_sa1 != npos);
    if (p_sb1 == npos || p_sb2 == npos || p_sa1 == npos) return;
    // qB's non-adjacent hits are merged into its one object, in input order.
    CHECK(p_qb < p_sb1);
    CHECK(p_sb1 < p_sb2);
    CHECK(p_sb2 < p_qa);
    CHECK(p_qa < p_sa1);
}

// A skip marker only turns into a "skipped" object when the query produced no
// hits at all; otherwise the query is reported as searched.
static void test_json_skip_and_hit_same_query() {
    std::fprintf(stderr, "-- test_json_skip_and_hit_same_query\n");

    OutputHit skip;
    skip.qseqid = "qA";
    skip.qlen = 30;
    skip.sstrand = '*';
    skip.skip_reason = ikafssn::kSkipQueryTooShort;
    skip.skip_detail = "too short";

    OutputHit hit_b;
    hit_b.qseqid = "qB"; hit_b.sseqid = "sB"; hit_b.sstrand = '+';
    hit_b.qstart = 0; hit_b.qend = 10; hit_b.qlen = 10;
    hit_b.sstart = 0; hit_b.send = 10; hit_b.slen = 100;
    hit_b.volume = 0;

    OutputHit hit_a = hit_b;
    hit_a.qseqid = "qA";
    hit_a.sseqid = "sA";

    const std::size_t npos = std::string::npos;
    {
        std::vector<OutputHit> hits{skip, hit_b, hit_a};
        std::ostringstream oss;
        write_results_json(oss, hits, /*mode=*/2, /*stage3_traceback=*/false);
        const std::string out = oss.str();
        const std::size_t p_qa = out.find("\"qseqid\": \"qA\"");
        const std::size_t p_qb = out.find("\"qseqid\": \"qB\"");
        CHECK(p_qa != npos);
        CHECK(p_qb != npos);
        if (p_qa == npos || p_qb == npos) return;
        // qA is first because its skip marker came first, and it reports the
        // hit it did produce.
        CHECK(p_qa < p_qb);
        CHECK(out.find("\"status\": \"skipped\"") == npos);
        CHECK(out.find("\"sseqid\": \"sA\"") != npos);
    }
    {
        std::vector<OutputHit> hits{skip, hit_b};
        std::ostringstream oss;
        write_results_json(oss, hits, /*mode=*/2, /*stage3_traceback=*/false);
        const std::string out = oss.str();
        const std::size_t p_qa = out.find("\"qseqid\": \"qA\"");
        const std::size_t p_qb = out.find("\"qseqid\": \"qB\"");
        const std::size_t p_sk = out.find("\"status\": \"skipped\"");
        CHECK(p_qa != npos);
        CHECK(p_qb != npos);
        CHECK(p_sk != npos);
        if (p_qa == npos || p_qb == npos || p_sk == npos) return;
        CHECK(p_qa < p_qb);
        CHECK(p_sk < p_qb);
    }
}

int main() {
    g_test_dir = "/tmp/ikafssn_result_reader_test";
    std::filesystem::create_directories(g_test_dir);

    test_basic_parse();
    test_skip_header_and_blank();
    test_invalid_lines();
    test_empty_input();
    test_header_only();
    test_stream_interface();
    test_roundtrip();
    test_roundtrip_mode3_no_traceback();
    test_roundtrip_mode3_traceback();
    test_numeric_field_formatting();
    test_json_string_escaping();
    test_json_query_first_appearance_order();
    test_json_skip_and_hit_same_query();
    test_writer_thread_count_invariance();
    test_writer_parallel_chunking();
    test_header_reordered_columns();
    test_no_header_fallback();
    test_has_volume_column();
    test_output_coords_roundtrip();
    test_result_reader_blast_coords();
    test_windows_line_endings();

    std::filesystem::remove_all(g_test_dir);

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
