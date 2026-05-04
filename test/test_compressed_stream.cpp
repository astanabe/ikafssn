// Round-trip + concatenation + magic detection tests for
// src/io/compressed_stream.{hpp,cpp}.

#include "io/compressed_stream.hpp"
#include "io/fasta_reader.hpp"
#include "io/result_writer.hpp"

#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "test_util.hpp"

using ikafssn::CompressionFormat;
using ikafssn::detect_format_from_extension;
using ikafssn::detect_format_from_magic;
using ikafssn::kCompressionLevelDefault;
using ikafssn::open_input_compressed;
using ikafssn::open_output_compressed;
using ikafssn::OutputFormat;
using ikafssn::validate_compression_level;
using ikafssn::validate_output_format;

namespace fs = std::filesystem;

static std::string test_dir() {
    auto dir = fs::temp_directory_path() / test_tmpdir("ikafssn_compressed_stream_test");
    std::error_code ec;
    fs::create_directories(dir, ec);
    return dir.string();
}

static std::string make_payload(std::size_t n, std::mt19937& rng) {
    static const char alphabet[] =
        ">ACGTNacgtn\n"
        "abcdefghijklmnopqrstuvwxyz0123456789 ";
    std::uniform_int_distribution<int> d(0, sizeof(alphabet) - 2);
    std::string s;
    s.reserve(n);
    for (std::size_t i = 0; i < n; i++) s.push_back(alphabet[d(rng)]);
    return s;
}

static bool write_via(const std::string& path, const std::string& data,
                      int level) {
    std::string err;
    auto out = open_output_compressed(path, level, err);
    if (!out) {
        std::fprintf(stderr, "open_output_compressed: %s\n", err.c_str());
        return false;
    }
    out.stream->write(data.data(), static_cast<std::streamsize>(data.size()));
    return out.stream->good();
}

static std::string read_via(const std::string& path) {
    std::string err;
    auto in = open_input_compressed(path, err);
    if (!in) {
        std::fprintf(stderr, "open_input_compressed: %s\n", err.c_str());
        return {};
    }
    std::ostringstream oss;
    oss << in.stream->rdbuf();
    return oss.str();
}

static void test_round_trip_each_codec() {
    std::mt19937 rng(0xC0DEC);
    auto dir = test_dir();
    struct Case { const char* name; const char* ext; };
    Case cases[] = {
        {"gzip",  "gz"},
        {"bzip2", "bz2"},
        {"xz",    "xz"},
        {"zstd",  "zst"},
    };
    for (auto& c : cases) {
        std::string payload = make_payload(64 * 1024 + 137, rng);
        std::string path = dir + "/round_trip_" + c.name + "." + c.ext;
        CHECK(write_via(path, payload, kCompressionLevelDefault));
        std::string back = read_via(path);
        CHECK(back == payload);
    }
}

static void test_concatenated_input() {
    std::mt19937 rng(0xC4F);
    auto dir = test_dir();
    struct Case { const char* name; const char* ext; };
    Case cases[] = {
        {"gzip",  "gz"},
        {"bzip2", "bz2"},
        {"xz",    "xz"},
        {"zstd",  "zst"},
    };
    for (auto& c : cases) {
        std::string a = make_payload(8 * 1024 + 11, rng);
        std::string b = make_payload(12 * 1024 + 23, rng);
        std::string p1 = dir + "/concat_a_" + c.name + "." + c.ext;
        std::string p2 = dir + "/concat_b_" + c.name + "." + c.ext;
        std::string p3 = dir + "/concat_ab_" + c.name + "." + c.ext;
        CHECK(write_via(p1, a, kCompressionLevelDefault));
        CHECK(write_via(p2, b, kCompressionLevelDefault));
        // cat p1 p2 > p3
        std::ifstream i1(p1, std::ios::binary);
        std::ifstream i2(p2, std::ios::binary);
        std::ofstream o3(p3, std::ios::binary);
        o3 << i1.rdbuf();
        o3 << i2.rdbuf();
        o3.close();
        std::string back = read_via(p3);
        CHECK(back == a + b);
    }
}

static void test_buffer_boundary_gzip() {
    // Pick a payload size large enough that the produced gzip stream
    // straddles the 64 KiB raw-input buffer of the decoder.  The exact
    // ratio is not important; ~256 KiB of randomised payload reliably
    // produces several 64 KiB blocks of compressed output.
    std::mt19937 rng(0xB0FF);
    std::string a = make_payload(256 * 1024 + 7, rng);
    std::string b = make_payload(192 * 1024 + 23, rng);

    auto dir = test_dir();
    std::string p1 = dir + "/buf_boundary_a.gz";
    std::string p2 = dir + "/buf_boundary_b.gz";
    std::string p3 = dir + "/buf_boundary_ab.gz";
    CHECK(write_via(p1, a, kCompressionLevelDefault));
    CHECK(write_via(p2, b, kCompressionLevelDefault));
    std::ifstream i1(p1, std::ios::binary);
    std::ifstream i2(p2, std::ios::binary);
    std::ofstream o3(p3, std::ios::binary);
    o3 << i1.rdbuf();
    o3 << i2.rdbuf();
    o3.close();
    std::string back = read_via(p3);
    CHECK(back == a + b);
}

static void test_magic_detection() {
    // Each codec's first few bytes should be classified correctly even
    // when the full prefix (6 bytes) is not present — except for xz,
    // which needs all 6 bytes of its magic.
    unsigned char gzip[]  = {0x1F, 0x8B};
    unsigned char bzip2[] = {'B', 'Z', 'h'};
    unsigned char xz[]    = {0xFD, '7', 'z', 'X', 'Z', 0x00};
    unsigned char zstd[]  = {0x28, 0xB5, 0x2F, 0xFD};
    unsigned char text[]  = {'>', 's', 'e', 'q', '1', '\n'};
    unsigned char short_prefix[] = {'>'};

    CHECK(detect_format_from_magic(gzip,  2) == CompressionFormat::kGzip);
    CHECK(detect_format_from_magic(bzip2, 3) == CompressionFormat::kBzip2);
    CHECK(detect_format_from_magic(xz,    6) == CompressionFormat::kXz);
    CHECK(detect_format_from_magic(zstd,  4) == CompressionFormat::kZstd);
    CHECK(detect_format_from_magic(text,  6) == CompressionFormat::kNone);
    CHECK(detect_format_from_magic(short_prefix, 1) == CompressionFormat::kNone);
}

static void test_extension_detection() {
    CHECK(detect_format_from_extension("foo.fa.gz")  == CompressionFormat::kGzip);
    CHECK(detect_format_from_extension("foo.tsv.bz2") == CompressionFormat::kBzip2);
    CHECK(detect_format_from_extension("foo.xz")     == CompressionFormat::kXz);
    CHECK(detect_format_from_extension("foo.zst")    == CompressionFormat::kZstd);
    // Case-insensitive
    CHECK(detect_format_from_extension("foo.GZ")     == CompressionFormat::kGzip);
    CHECK(detect_format_from_extension("foo.ZsT")    == CompressionFormat::kZstd);
    // Plain / unrelated suffixes
    CHECK(detect_format_from_extension("foo.fa")     == CompressionFormat::kNone);
    CHECK(detect_format_from_extension("foo")        == CompressionFormat::kNone);
    CHECK(detect_format_from_extension("")           == CompressionFormat::kNone);
    CHECK(detect_format_from_extension("-")          == CompressionFormat::kNone);
}

static void test_level_validation() {
    std::string err;
    // Default sentinel always accepted
    CHECK(validate_compression_level(CompressionFormat::kGzip,
                                      kCompressionLevelDefault, err));
    CHECK(validate_compression_level(CompressionFormat::kBzip2,
                                      kCompressionLevelDefault, err));
    CHECK(validate_compression_level(CompressionFormat::kXz,
                                      kCompressionLevelDefault, err));
    CHECK(validate_compression_level(CompressionFormat::kZstd,
                                      kCompressionLevelDefault, err));

    // gzip 0..9
    CHECK(validate_compression_level(CompressionFormat::kGzip, 0, err));
    CHECK(validate_compression_level(CompressionFormat::kGzip, 9, err));
    err.clear();
    CHECK(!validate_compression_level(CompressionFormat::kGzip, 10, err));
    CHECK(err.find("gzip") != std::string::npos);
    CHECK(err.find("0..9") != std::string::npos);

    // bzip2 1..9 (rejects 0)
    err.clear();
    CHECK(!validate_compression_level(CompressionFormat::kBzip2, 0, err));
    CHECK(validate_compression_level(CompressionFormat::kBzip2, 1, err));
    CHECK(validate_compression_level(CompressionFormat::kBzip2, 9, err));

    // xz 0..9
    CHECK(validate_compression_level(CompressionFormat::kXz, 0, err));
    CHECK(validate_compression_level(CompressionFormat::kXz, 9, err));

    // zstd has a wide range; just confirm 1 is acceptable and that some
    // wildly-out value is rejected.
    err.clear();
    CHECK(validate_compression_level(CompressionFormat::kZstd, 1, err));
    err.clear();
    CHECK(!validate_compression_level(CompressionFormat::kZstd, 1000000, err));
    CHECK(err.find("zstd") != std::string::npos);
}

static void test_sam_bam_compression_rejection() {
    std::string err;
    // SAM with .gz suffix → rejected
    err.clear();
    CHECK(!validate_output_format(OutputFormat::kSam, /*mode=*/3,
                                   /*traceback=*/true, "out.sam.gz", err));
    CHECK(err.find("compression") != std::string::npos);

    err.clear();
    CHECK(!validate_output_format(OutputFormat::kBam, /*mode=*/3,
                                   /*traceback=*/true, "out.bam.zst", err));

    // Plain SAM still allowed
    err.clear();
    CHECK(validate_output_format(OutputFormat::kSam, /*mode=*/3,
                                  /*traceback=*/true, "out.sam", err));
    // Plain BAM still allowed
    err.clear();
    CHECK(validate_output_format(OutputFormat::kBam, /*mode=*/3,
                                  /*traceback=*/true, "out.bam", err));

    // TSV with .gz is fine
    err.clear();
    CHECK(validate_output_format(OutputFormat::kTab, /*mode=*/2,
                                  /*traceback=*/false, "out.tsv.gz", err));
}

static void test_fasta_round_trip_through_gzip() {
    auto dir = test_dir();
    std::string plain = dir + "/queries.fa";
    std::string gz    = dir + "/queries.fa.gz";

    std::string body =
        ">seqA\nACGTACGTACGT\n"
        ">seqB\nGGGGCCCCAAAATTTT\n"
        ">seqC\nNNNNACGTNNNN\n";

    // Plain
    {
        std::ofstream o(plain);
        o << body;
    }
    // gzipped
    CHECK(write_via(gz, body, kCompressionLevelDefault));

    auto plain_recs = ikafssn::read_fasta(plain);
    auto gz_recs    = ikafssn::read_fasta(gz);
    CHECK(plain_recs.size() == 3);
    CHECK(gz_recs.size() == 3);
    CHECK(plain_recs.size() == gz_recs.size());
    for (std::size_t i = 0; i < plain_recs.size(); i++) {
        CHECK(plain_recs[i].id == gz_recs[i].id);
        CHECK(plain_recs[i].sequence == gz_recs[i].sequence);
    }
}

int main() {
    test_round_trip_each_codec();
    test_concatenated_input();
    test_buffer_boundary_gzip();
    test_magic_detection();
    test_extension_detection();
    test_level_validation();
    test_sam_bam_compression_rejection();
    test_fasta_round_trip_through_gzip();
    TEST_SUMMARY();
    return g_fail_count == 0 ? 0 : 1;
}
