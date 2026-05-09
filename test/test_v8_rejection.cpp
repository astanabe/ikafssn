#include "test_util.hpp"
#include "index/kix_reader.hpp"
#include "index/kpx_reader.hpp"
#include "index/ksx_reader.hpp"
#include "index/khx_reader.hpp"
#include "index/kix_format.hpp"
#include "index/kpx_format.hpp"
#include "index/ksx_format.hpp"
#include "index/khx_format.hpp"
#include "core/config.hpp"

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

using namespace ikafssn;

namespace {

// Synthesize a minimal .kix file with v8 magic + format_version 8 and
// verify KixReader::open rejects it.  We do not have to build a fully
// valid v8 body — KixReader rejects on the magic / version check
// before reading anything else.
std::string synth_v8_kix(const std::string& path) {
    KixHeader hdr{};
    // v8 magic
    hdr.magic[0] = 'K'; hdr.magic[1] = 'I'; hdr.magic[2] = 'X'; hdr.magic[3] = '8';
    hdr.format_version = 8;
    hdr.k = 5;
    FILE* fp = std::fopen(path.c_str(), "wb");
    std::fwrite(&hdr, sizeof(hdr), 1, fp);
    // Pad with a few zero bytes so the file is at least header-sized.
    std::vector<uint8_t> filler(128, 0);
    std::fwrite(filler.data(), 1, filler.size(), fp);
    std::fclose(fp);
    return path;
}

std::string synth_v8_kpx(const std::string& path) {
    KpxHeader hdr{};
    hdr.magic[0] = 'K'; hdr.magic[1] = 'P'; hdr.magic[2] = 'X'; hdr.magic[3] = '8';
    hdr.format_version = 8;
    hdr.k = 5;
    FILE* fp = std::fopen(path.c_str(), "wb");
    std::fwrite(&hdr, sizeof(hdr), 1, fp);
    std::vector<uint8_t> filler(128, 0);
    std::fwrite(filler.data(), 1, filler.size(), fp);
    std::fclose(fp);
    return path;
}

std::string synth_v8_ksx(const std::string& path) {
    KsxHeader hdr{};
    std::memcpy(hdr.magic, KSX_MAGIC, 4);  // KMSX (unchanged)
    hdr.format_version = 8;
    FILE* fp = std::fopen(path.c_str(), "wb");
    std::fwrite(&hdr, sizeof(hdr), 1, fp);
    std::fclose(fp);
    return path;
}

std::string synth_v8_khx(const std::string& path) {
    KhxHeader hdr{};
    std::memcpy(hdr.magic, KHX_MAGIC, 4);  // KMHX (unchanged)
    hdr.format_version = 8;
    hdr.k = 5;
    FILE* fp = std::fopen(path.c_str(), "wb");
    std::fwrite(&hdr, sizeof(hdr), 1, fp);
    std::fclose(fp);
    return path;
}

void test_kix_v8_rejected() {
    std::string path = test_tmpdir("/tmp/test_v8_kix") + ".kix";
    synth_v8_kix(path);
    KixReader reader;
    CHECK(!reader.open(path));
    std::remove(path.c_str());
}

void test_kpx_v8_rejected() {
    std::string path = test_tmpdir("/tmp/test_v8_kpx") + ".kpx";
    synth_v8_kpx(path);
    KpxReader reader;
    CHECK(!reader.open(path));
    std::remove(path.c_str());
}

void test_ksx_v8_rejected() {
    std::string path = test_tmpdir("/tmp/test_v8_ksx") + ".ksx";
    synth_v8_ksx(path);
    KsxReader reader;
    CHECK(!reader.open(path));
    std::remove(path.c_str());
}

void test_khx_v8_rejected() {
    std::string path = test_tmpdir("/tmp/test_v8_khx") + ".khx";
    synth_v8_khx(path);
    KhxReader reader;
    CHECK(!reader.open(path));
    std::remove(path.c_str());
}

// Confirm v10 is the *current* format constant.  Guards against an
// accidental constant downgrade.
void test_v10_is_current() {
    CHECK_EQ(KIX_FORMAT_VERSION, 10u);
    CHECK_EQ(KPX_FORMAT_VERSION, 10u);
    CHECK_EQ(KSX_FORMAT_VERSION, 10u);
    CHECK_EQ(KHX_FORMAT_VERSION, 10u);
    // v10 magic widened from 4 bytes "KIX9" to 5 bytes "KIX10".
    CHECK(KIX_MAGIC[3] == '1' && KIX_MAGIC[4] == '0');
    CHECK(KPX_MAGIC[3] == '1' && KPX_MAGIC[4] == '0');
}

// Synthesize a .kix file with the *wrong* magic ("KIX08" — the v10
// 5-byte magic field but with v9-style version digit) while format
// version is 10.  The reader should reject on magic mismatch first.
void test_wrong_magic_rejected() {
    std::string path = test_tmpdir("/tmp/test_v8_wrong_magic") + ".kix";
    KixHeader hdr{};
    hdr.magic[0] = 'K'; hdr.magic[1] = 'I'; hdr.magic[2] = 'X';
    hdr.magic[3] = '0'; hdr.magic[4] = '8';
    hdr.format_version = KIX_FORMAT_VERSION;  // = 10
    hdr.k = 5;
    FILE* fp = std::fopen(path.c_str(), "wb");
    std::fwrite(&hdr, sizeof(hdr), 1, fp);
    std::vector<uint8_t> filler(128, 0);
    std::fwrite(filler.data(), 1, filler.size(), fp);
    std::fclose(fp);
    KixReader reader;
    CHECK(!reader.open(path));
    std::remove(path.c_str());
}

} // anonymous namespace

int main() {
    test_kix_v8_rejected();
    test_kpx_v8_rejected();
    test_ksx_v8_rejected();
    test_khx_v8_rejected();
    test_v10_is_current();
    test_wrong_magic_rejected();
    TEST_SUMMARY();
    return g_fail_count == 0 ? 0 : 1;
}
