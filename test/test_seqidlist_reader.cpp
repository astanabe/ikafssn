#include "test_util.hpp"
#include "io/seqidlist_reader.hpp"

#include <fstream>
#include <string>
#include <filesystem>

using namespace ikafssn;

static std::string g_test_dir;

static void write_file(const std::string& path, const std::string& content) {
    std::ofstream f(path);
    f << content;
}

static void test_basic_text() {
    std::fprintf(stderr, "-- test_basic_text\n");

    std::string path = g_test_dir + "/basic.txt";
    write_file(path, "ACC001\nACC002\nACC003\n");

    auto result = read_seqidlist(path);
    CHECK_EQ(result.size(), 3u);
    CHECK(result[0] == "ACC001");
    CHECK(result[1] == "ACC002");
    CHECK(result[2] == "ACC003");
}

static void test_trim_angle_bracket() {
    std::fprintf(stderr, "-- test_trim_angle_bracket\n");

    std::string path = g_test_dir + "/bracket.txt";
    write_file(path, ">ACC001\n>ACC002\n");

    auto result = read_seqidlist(path);
    CHECK_EQ(result.size(), 2u);
    CHECK(result[0] == "ACC001");
    CHECK(result[1] == "ACC002");
}

static void test_skip_blank_lines() {
    std::fprintf(stderr, "-- test_skip_blank_lines\n");

    std::string path = g_test_dir + "/blank.txt";
    write_file(path, "ACC001\n\n   \n\nACC002\n");

    auto result = read_seqidlist(path);
    CHECK_EQ(result.size(), 2u);
}

static void test_skip_comments() {
    std::fprintf(stderr, "-- test_skip_comments\n");

    std::string path = g_test_dir + "/comments.txt";
    write_file(path, "# this is a comment\nACC001\n# another comment\nACC002\n");

    auto result = read_seqidlist(path);
    CHECK_EQ(result.size(), 2u);
}

static void test_first_token_only() {
    std::fprintf(stderr, "-- test_first_token_only\n");

    std::string path = g_test_dir + "/tokens.txt";
    write_file(path, "ACC001 some extra text\nACC002\tmore\n");

    auto result = read_seqidlist(path);
    CHECK_EQ(result.size(), 2u);
    CHECK(result[0] == "ACC001");
    CHECK(result[1] == "ACC002");
}

static void test_empty_file() {
    std::fprintf(stderr, "-- test_empty_file\n");

    std::string path = g_test_dir + "/empty.txt";
    write_file(path, "");

    auto result = read_seqidlist(path);
    CHECK_EQ(result.size(), 0u);
}

static void test_nonexistent_file() {
    std::fprintf(stderr, "-- test_nonexistent_file\n");

    auto result = read_seqidlist("/nonexistent/path/file.txt");
    CHECK_EQ(result.size(), 0u);
}

// --- Binary seqidlist format tests ---

// Helper to write a little-endian uint64
static void write_u64(std::ofstream& f, uint64_t v) {
    unsigned char buf[8];
    for (int i = 0; i < 8; i++) { buf[i] = v & 0xFF; v >>= 8; }
    f.write(reinterpret_cast<char*>(buf), 8);
}

// Helper to write a little-endian uint32
static void write_u32(std::ofstream& f, uint32_t v) {
    unsigned char buf[4];
    for (int i = 0; i < 4; i++) { buf[i] = v & 0xFF; v >>= 8; }
    f.write(reinterpret_cast<char*>(buf), 4);
}

// Create a minimal v5 binary seqidlist file
static void create_binary_seqidlist(const std::string& path,
                                    const std::vector<std::string>& ids) {
    // Calculate data section size
    size_t data_size = 0;
    for (const auto& id : ids) {
        if (id.size() >= 255) {
            data_size += 1 + 4 + id.size(); // 0xFF marker + uint32 len + string
        } else {
            data_size += 1 + id.size(); // uint8 len + string
        }
    }

    std::string title = "test";
    std::string create_date = "2026-01-01";
    // Header size: 1 (null) + 8 (file_size) + 8 (num_ids) + 4 (title_len) + title
    //            + 1 (date_len) + date + 8 (db_vol_length=0)
    size_t header_size = 1 + 8 + 8 + 4 + title.size() + 1 + create_date.size() + 8;
    size_t total_size = header_size + data_size;

    std::ofstream f(path, std::ios::binary);

    // Null byte (v5 marker)
    char zero = 0x00;
    f.write(&zero, 1);

    // File size
    write_u64(f, total_size);

    // Number of IDs
    write_u64(f, ids.size());

    // Title
    write_u32(f, static_cast<uint32_t>(title.size()));
    f.write(title.data(), title.size());

    // Create date
    unsigned char date_len = static_cast<unsigned char>(create_date.size());
    f.write(reinterpret_cast<char*>(&date_len), 1);
    f.write(create_date.data(), create_date.size());

    // DB volume length = 0 (not associated with DB)
    write_u64(f, 0);

    // Data section: IDs
    for (const auto& id : ids) {
        if (id.size() >= 255) {
            unsigned char marker = 0xFF;
            f.write(reinterpret_cast<char*>(&marker), 1);
            write_u32(f, static_cast<uint32_t>(id.size()));
        } else {
            unsigned char len = static_cast<unsigned char>(id.size());
            f.write(reinterpret_cast<char*>(&len), 1);
        }
        f.write(id.data(), id.size());
    }
}

static void test_binary_basic() {
    std::fprintf(stderr, "-- test_binary_basic\n");

    std::string path = g_test_dir + "/basic.bsl";
    create_binary_seqidlist(path, {"NM_001234", "XM_005678", "NR_999999"});

    auto result = read_seqidlist(path);
    CHECK_EQ(result.size(), 3u);
    CHECK(result[0] == "NM_001234");
    CHECK(result[1] == "XM_005678");
    CHECK(result[2] == "NR_999999");
}

static void test_binary_single_id() {
    std::fprintf(stderr, "-- test_binary_single_id\n");

    std::string path = g_test_dir + "/single.bsl";
    create_binary_seqidlist(path, {"ACC001"});

    auto result = read_seqidlist(path);
    CHECK_EQ(result.size(), 1u);
    CHECK(result[0] == "ACC001");
}

static void test_binary_empty() {
    std::fprintf(stderr, "-- test_binary_empty\n");

    std::string path = g_test_dir + "/empty.bsl";
    create_binary_seqidlist(path, {});

    auto result = read_seqidlist(path);
    CHECK_EQ(result.size(), 0u);
}

static void test_binary_long_id() {
    std::fprintf(stderr, "-- test_binary_long_id\n");

    // Create an ID longer than 255 bytes to test the 0xFF + uint32 path
    std::string long_id(300, 'A');
    std::string path = g_test_dir + "/long.bsl";
    create_binary_seqidlist(path, {"SHORT", long_id, "AFTER"});

    auto result = read_seqidlist(path);
    CHECK_EQ(result.size(), 3u);
    CHECK(result[0] == "SHORT");
    CHECK(result[1] == long_id);
    CHECK(result[2] == "AFTER");
}

static void test_binary_autodetect() {
    std::fprintf(stderr, "-- test_binary_autodetect\n");

    // Verify text files are NOT misdetected as binary
    std::string text_path = g_test_dir + "/not_binary.txt";
    write_file(text_path, "ACC001\nACC002\n");
    auto text_result = read_seqidlist(text_path);
    CHECK_EQ(text_result.size(), 2u);

    // Verify binary files ARE detected as binary
    std::string bin_path = g_test_dir + "/is_binary.bsl";
    create_binary_seqidlist(bin_path, {"ACC001", "ACC002"});
    auto bin_result = read_seqidlist(bin_path);
    CHECK_EQ(bin_result.size(), 2u);

    // Both should produce the same accessions
    CHECK(text_result[0] == bin_result[0]);
    CHECK(text_result[1] == bin_result[1]);
}

int main() {
    g_test_dir = "/tmp/ikafssn_seqidlist_test";
    std::filesystem::create_directories(g_test_dir);

    // Text format tests
    test_basic_text();
    test_trim_angle_bracket();
    test_skip_blank_lines();
    test_skip_comments();
    test_first_token_only();
    test_empty_file();
    test_nonexistent_file();

    // Binary format tests
    test_binary_basic();
    test_binary_single_id();
    test_binary_empty();
    test_binary_long_id();
    test_binary_autodetect();

    std::filesystem::remove_all(g_test_dir);

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
