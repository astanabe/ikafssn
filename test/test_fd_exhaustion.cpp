// ikafssnindex and ikafssnretrieve under a tight RLIMIT_NOFILE.
//
// Both walk the BLAST DB one volume at a time, so their descriptor usage is
// a small constant no matter how many volumes the database has.  Each is run
// in a child process whose descriptor limit is set well below what opening
// every volume at once would take.

#include "test_util.hpp"
#include "ssu_test_fixture.hpp"
#include "io/blastdb_reader.hpp"
#include "io/result_writer.hpp"

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>

using namespace ikafssn;
using namespace ssu_fixture;

// Descriptor budget handed to the children.  Opening every volume of
// ssu_manyvol at once costs three descriptors per volume (.nin + .nsq +
// .nhr), which overruns this; opening them one at a time costs a handful.
static constexpr rlim_t kChildFdLimit = 40;

static std::string g_db_prefix;
static std::string g_test_dir;
static size_t g_num_volumes = 0;

// Run argv in a child capped at `kChildFdLimit` descriptors and return its
// exit status (-1 if it did not exit normally).
static int run_capped(const std::vector<std::string>& argv) {
    pid_t pid = ::fork();
    if (pid < 0) return -1;

    if (pid == 0) {
        struct rlimit rl;
        rl.rlim_cur = kChildFdLimit;
        rl.rlim_max = kChildFdLimit;
        if (::setrlimit(RLIMIT_NOFILE, &rl) != 0) ::_exit(127);

        std::vector<char*> args;
        args.reserve(argv.size() + 1);
        for (const auto& a : argv) args.push_back(const_cast<char*>(a.c_str()));
        args.push_back(nullptr);
        ::execv(args[0], args.data());
        ::_exit(127);
    }

    int status = 0;
    if (::waitpid(pid, &status, 0) < 0) return -1;
    return WIFEXITED(status) ? WEXITSTATUS(status) : -1;
}

static void test_index_under_fd_limit() {
    std::fprintf(stderr, "-- test_index_under_fd_limit (%zu volumes, limit %llu)\n",
                 g_num_volumes, static_cast<unsigned long long>(kChildFdLimit));

    const std::string out_dir = g_test_dir + "/index";
    std::filesystem::create_directories(out_dir);

    int rc = run_capped({IKAFSSN_INDEX_BIN, "-db", g_db_prefix,
                         "-k", "7", "-o", out_dir});
    CHECK_EQ(rc, 0);

    // One .ksx per volume proves every volume was indexed, not just the ones
    // that fitted under the limit.
    size_t ksx_count = 0;
    for (const auto& e : std::filesystem::directory_iterator(out_dir)) {
        if (e.path().extension() == ".ksx") ksx_count++;
    }
    CHECK_EQ(ksx_count, g_num_volumes);
}

// Write one hit per volume, taking the first OID of each, so retrieval has to
// visit every volume of the database.  The rows go out in descending volume
// order — the opposite of the order the volumes are walked in — so the output
// can only match `expected` if the records were staged and put back in input
// order.
static bool write_hits_tsv(const std::string& path,
                           std::vector<std::string>& expected) {
    auto vol_paths = BlastDbReader::find_volume_paths(g_db_prefix);
    std::vector<OutputHit> hits;

    for (size_t vi = vol_paths.size(); vi-- > 0; ) {
        BlastDbReader reader;
        if (!reader.open(vol_paths[vi])) return false;
        if (reader.num_sequences() == 0) continue;

        std::string acc;
        if (!reader.get_accession(0, acc)) return false;
        const uint32_t slen = reader.seq_length(0);
        if (acc.empty() || slen == 0) return false;

        OutputHit h;
        h.qseqid = "fdtest";
        h.sseqid = acc;
        h.sstrand = '+';
        h.qstart = 0;
        h.qend = slen;
        h.qlen = slen;
        h.sstart = 0;
        h.send = slen;
        h.slen = slen;
        h.chainscore = 1;
        h.volume = static_cast<uint16_t>(vi);
        expected.push_back(acc);
        hits.push_back(std::move(h));
    }

    if (hits.size() != vol_paths.size()) return false;

    std::ofstream ofs(path);
    if (!ofs) return false;
    write_results_tsv(ofs, hits);
    return static_cast<bool>(ofs);
}

// Copy `src` to `dst` with the last tab-separated field of every line removed.
// `volume` is the last column of every TSV layout ikafssn writes.
static bool strip_last_column(const std::string& src, const std::string& dst) {
    std::ifstream in(src);
    std::ofstream out(dst);
    if (!in || !out) return false;
    std::string line;
    while (std::getline(in, line)) {
        auto tab = line.rfind('\t');
        if (tab == std::string::npos) return false;
        out << line.substr(0, tab) << '\n';
    }
    return static_cast<bool>(out);
}

// Leading token of each FASTA defline, with the `:start-end` suffix removed.
static std::vector<std::string> defline_accessions(const std::string& fasta) {
    std::vector<std::string> out;
    std::ifstream ifs(fasta);
    std::string line;
    while (std::getline(ifs, line)) {
        if (line.empty() || line[0] != '>') continue;
        std::string tok = line.substr(1, line.find(' ') - 1);
        auto colon = tok.rfind(':');
        if (colon != std::string::npos) tok.resize(colon);
        out.push_back(std::move(tok));
    }
    return out;
}

static void test_retrieve_under_fd_limit() {
    std::fprintf(stderr, "-- test_retrieve_under_fd_limit\n");

    const std::string tsv = g_test_dir + "/hits.tsv";
    const std::string fasta = g_test_dir + "/hits.fasta";
    std::vector<std::string> expected;
    if (!write_hits_tsv(tsv, expected)) {
        std::fprintf(stderr, "  cannot build the hit list\n");
        CHECK(false);
        return;
    }

    int rc = run_capped({IKAFSSN_RETRIEVE_BIN, "-db", g_db_prefix,
                         "-tsv", tsv, "-o", fasta, "-context_extend", "0"});
    CHECK_EQ(rc, 0);

    // One record per volume, in input order rather than volume order.
    auto got = defline_accessions(fasta);
    CHECK_EQ(got.size(), g_num_volumes);
    CHECK(got == expected);

    // The scratch file is unlinked as soon as it is created, so the output
    // directory holds nothing but the inputs and the FASTA.
    for (const auto& e : std::filesystem::directory_iterator(g_test_dir)) {
        CHECK(e.path().filename().string().rfind("ikafssnretrieve.", 0) != 0);
    }
}

// Without a volume column there is no way to tell which BLAST DB volume a hit
// belongs to, so -db retrieval refuses the file rather than sending every hit
// to volume 0.
static void test_retrieve_rejects_missing_volume_column() {
    std::fprintf(stderr, "-- test_retrieve_rejects_missing_volume_column\n");

    const std::string tsv = g_test_dir + "/hits.tsv";
    const std::string stripped = g_test_dir + "/hits_novolume.tsv";
    const std::string fasta = g_test_dir + "/novolume.fasta";
    if (!strip_last_column(tsv, stripped)) {
        std::fprintf(stderr, "  cannot strip the volume column\n");
        CHECK(false);
        return;
    }

    int rc = run_capped({IKAFSSN_RETRIEVE_BIN, "-db", g_db_prefix,
                         "-tsv", stripped, "-o", fasta});
    CHECK_EQ(rc, 1);
    CHECK(!std::filesystem::exists(fasta) ||
          std::filesystem::file_size(fasta) == 0);
}

int main() {
    g_db_prefix = manyvol_prefix();
    if (!std::filesystem::exists(g_db_prefix + ".nal")) {
        skip("ssu_manyvol not found; run test/scripts/setup_ssu_testdata.sh");
    }

    g_num_volumes = BlastDbReader::find_volume_paths(g_db_prefix).size();
    if (g_num_volumes * 3 <= kChildFdLimit) {
        skip("ssu_manyvol has too few volumes to exercise the limit");
    }

    g_test_dir = test_tmpdir("/tmp/ikafssn_fd_exhaustion_test");
    std::filesystem::remove_all(g_test_dir);
    std::filesystem::create_directories(g_test_dir);

    test_index_under_fd_limit();
    test_retrieve_under_fd_limit();
    test_retrieve_rejects_missing_volume_column();

    std::filesystem::remove_all(g_test_dir);

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
