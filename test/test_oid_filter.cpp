#include "test_util.hpp"
#include "search/oid_filter.hpp"
#include "index/ksx_reader.hpp"
#include "index/ksx_writer.hpp"
#include "index/ksx_format.hpp"

#include <string>
#include <vector>
#include <filesystem>

using namespace ikafssn;

static std::string g_test_dir;

// Create a minimal .ksx file with known accessions
static void create_test_ksx(const std::string& path,
                            const std::vector<std::string>& accessions,
                            const std::vector<uint32_t>& lengths) {
    KsxWriter writer;
    for (size_t i = 0; i < accessions.size(); i++) {
        writer.add_sequence(lengths[i], accessions[i]);
    }
    writer.write(path);
}

static void test_no_filter() {
    std::fprintf(stderr, "-- test_no_filter\n");

    OidFilter filter;
    // Default mode is kNone, all should pass
    CHECK(filter.pass(0));
    CHECK(filter.pass(1));
    CHECK(filter.pass(100));
}

static void test_include_mode() {
    std::fprintf(stderr, "-- test_include_mode\n");

    std::string ksx_path = g_test_dir + "/include.ksx";
    create_test_ksx(ksx_path, {"ACC1", "ACC2", "ACC3", "ACC4"}, {100, 200, 300, 400});

    KsxReader ksx;
    CHECK(ksx.open(ksx_path));

    OidFilter filter;
    filter.build({"ACC1", "ACC3"}, ksx, OidFilterMode::kInclude);

    CHECK(filter.pass(0));   // ACC1 - included
    CHECK(!filter.pass(1));  // ACC2 - not in list
    CHECK(filter.pass(2));   // ACC3 - included
    CHECK(!filter.pass(3));  // ACC4 - not in list

    ksx.close();
}

static void test_exclude_mode() {
    std::fprintf(stderr, "-- test_exclude_mode\n");

    std::string ksx_path = g_test_dir + "/exclude.ksx";
    create_test_ksx(ksx_path, {"ACC1", "ACC2", "ACC3", "ACC4"}, {100, 200, 300, 400});

    KsxReader ksx;
    CHECK(ksx.open(ksx_path));

    OidFilter filter;
    filter.build({"ACC2", "ACC4"}, ksx, OidFilterMode::kExclude);

    CHECK(filter.pass(0));   // ACC1 - not excluded
    CHECK(!filter.pass(1));  // ACC2 - excluded
    CHECK(filter.pass(2));   // ACC3 - not excluded
    CHECK(!filter.pass(3));  // ACC4 - excluded

    ksx.close();
}

static void test_unresolved_accession() {
    std::fprintf(stderr, "-- test_unresolved_accession\n");

    std::string ksx_path = g_test_dir + "/unresolved.ksx";
    create_test_ksx(ksx_path, {"ACC1", "ACC2"}, {100, 200});

    KsxReader ksx;
    CHECK(ksx.open(ksx_path));

    OidFilter filter;
    // "UNKNOWN" won't be found; should warn but not crash
    filter.build({"ACC1", "UNKNOWN"}, ksx, OidFilterMode::kInclude);

    CHECK(filter.pass(0));   // ACC1 - included
    CHECK(!filter.pass(1));  // ACC2 - not in list

    ksx.close();
}

static void test_empty_accessions() {
    std::fprintf(stderr, "-- test_empty_accessions\n");

    std::string ksx_path = g_test_dir + "/empty_acc.ksx";
    create_test_ksx(ksx_path, {"ACC1", "ACC2"}, {100, 200});

    KsxReader ksx;
    CHECK(ksx.open(ksx_path));

    OidFilter filter;
    filter.build({}, ksx, OidFilterMode::kInclude);

    // Empty accession list with include mode → reverts to kNone
    CHECK(filter.pass(0));
    CHECK(filter.pass(1));

    ksx.close();
}

int main() {
    g_test_dir = "/tmp/ikafssn_oid_filter_test";
    std::filesystem::create_directories(g_test_dir);

    test_no_filter();
    test_include_mode();
    test_exclude_mode();
    test_unresolved_accession();
    test_empty_accessions();

    std::filesystem::remove_all(g_test_dir);

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
