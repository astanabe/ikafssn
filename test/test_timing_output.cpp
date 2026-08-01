// bench/run_warm_e2e.sh attributes a run by parsing the key=value pairs of
// both timing lines ikafssnsearch prints at info level, so pin their key sets.
// Renaming one key would silently drop a metric from every measurement instead
// of failing loudly.  The "Timing run_search" line comes from the shared
// orchestrator and is pinned at the library level by test_multivolume; only the
// CLI emits "Timing overall", so it takes a run of the real binary.
#include "test_util.hpp"
#include "ssu_test_fixture.hpp"
#include "io/blastdb_reader.hpp"

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

using namespace ikafssn;
using namespace ssu_fixture;

namespace {

std::string g_build_dir;
std::string g_work_dir;
std::string g_index_prefix;
std::string g_query_path;

// Run a command with stderr merged into stdout and return everything it wrote.
std::string run_capture(const std::string& cmd) {
    FILE* pipe = popen((cmd + " 2>&1").c_str(), "r");
    if (!pipe) return {};
    std::string out;
    char buf[4096];
    while (std::fgets(buf, sizeof(buf), pipe)) out += buf;
    pclose(pipe);
    return out;
}

// Return the last line of `text` containing `needle`, or an empty string.
std::string last_line_with(const std::string& text, const std::string& needle) {
    std::istringstream in(text);
    std::string line, found;
    while (std::getline(in, line)) {
        if (line.find(needle) != std::string::npos) found = line;
    }
    return found;
}

void setup() {
    check_ssu_available();

    g_build_dir = std::string(SOURCE_DIR) + "/build/src";
    if (!std::filesystem::exists(g_build_dir + "/ikafssnindex") ||
        !std::filesystem::exists(g_build_dir + "/ikafssnsearch")) {
        skip("ikafssnindex / ikafssnsearch binaries not found (expected under build/src/)");
    }

    g_work_dir = test_tmpdir("/tmp/ikafssn_timing_test");
    std::error_code ec;
    std::filesystem::remove_all(g_work_dir, ec);
    std::filesystem::create_directories(g_work_dir);

    const std::string db = ssu_db_prefix();

    // A mode 1 index is enough: the "Timing overall" key set does not depend
    // on the mode.
    const std::string build_cmd =
        g_build_dir + "/ikafssnindex -db " + db + " -k 11 -o " + g_work_dir +
        " -mode 1 -max_degen_expand 0";
    const std::string build_log = run_capture(build_cmd);
    CHECK(build_log.find("All volumes completed successfully.") != std::string::npos);

    g_index_prefix = g_work_dir + "/" +
                     std::filesystem::path(db).filename().string();

    // Take the query out of the DB so the test needs no derived fixture data.
    BlastDbReader reader;
    CHECK(reader.open(db));
    const std::string seq = extract_subsequence(reader, ACC_FJ, 0, 300);
    CHECK(!seq.empty());

    g_query_path = g_work_dir + "/query.fasta";
    std::ofstream q(g_query_path);
    q << ">" << ACC_FJ << "\n" << seq << "\n";
    q.close();
}

void test_timing_line_keys() {
    std::fprintf(stderr, "-- test_timing_line_keys\n");

    const std::string out_path = g_work_dir + "/out.tsv";
    const std::string log = run_capture(
        g_build_dir + "/ikafssnsearch -ix " + g_index_prefix +
        " -k 11 -query " + g_query_path + " -o " + out_path + " -mode 1");

    const std::string overall = last_line_with(log, "Timing overall (s):");
    CHECK(!overall.empty());
    for (const char* key : {"open_index=", "preprocess=", "run_search=",
                            "convert=", "stage3=", "stage3_select=", "write=",
                            "total="}) {
        CHECK(overall.find(key) != std::string::npos);
    }

    // The CLI reaches the orchestrator's line too, so the two stay in sync.
    CHECK(!last_line_with(log, "Timing run_search (s):").empty());
}

void cleanup() {
    std::error_code ec;
    std::filesystem::remove_all(g_work_dir, ec);
}

} // namespace

int main() {
    setup();
    test_timing_line_keys();
    cleanup();
    TEST_SUMMARY();
    return g_fail_count == 0 ? 0 : 1;
}
