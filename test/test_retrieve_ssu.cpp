//
// Integration test for ikafssnretrieve using SSU_eukaryote_rRNA BLAST DB
// from NCBI FTP.
//
// This test downloads the DB (if not already present), then:
//   1. Tests local BLAST DB extraction for 3 known accessions
//   2. Tests NCBI efetch remote retrieval for the same accessions
//   3. Compares local and remote results for identity
//
// Requires: network access, NCBI Toolkit (local), libcurl (remote)
//

#include "test_util.hpp"
#include "io/blastdb_reader.hpp"
#include "io/result_writer.hpp"
#include "ikafssnretrieve/local_retriever.hpp"

#ifdef IKAFSSN_ENABLE_REMOTE
#include "ikafssnretrieve/efetch_retriever.hpp"
#endif

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace ikafssn;

static const char* SSU_URL = "https://ftp.ncbi.nih.gov/blast/db/SSU_eukaryote_rRNA.tar.gz";

// Target accessions (unversioned, as stored in BLAST DB)
static const char* TARGET_ACC[] = {"FJ876973", "GQ912721", "DQ235612"};
static constexpr int NUM_TARGETS = 3;

static std::string g_db_dir;
static std::string g_db_prefix;

// Per-accession info discovered from the BLAST DB
struct AccInfo {
    uint32_t seq_length = 0;
    bool found = false;
};
static AccInfo g_acc_info[NUM_TARGETS];

// ---- Helpers ----

// Parse FASTA from a string into accession -> sequence map.
// Header format: >ACCESSION query=... strand=... range=... score=...
static std::unordered_map<std::string, std::string>
parse_fasta_output(const std::string& fasta_str) {
    std::unordered_map<std::string, std::string> seqs;
    std::istringstream iss(fasta_str);
    std::string line, cur_acc, cur_seq;

    while (std::getline(iss, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!cur_acc.empty()) seqs[cur_acc] = cur_seq;
            // Extract first word after '>'
            size_t s = 1;
            while (s < line.size() && std::isspace(static_cast<unsigned char>(line[s]))) s++;
            size_t e = s;
            while (e < line.size() && !std::isspace(static_cast<unsigned char>(line[e]))) e++;
            cur_acc = line.substr(s, e - s);
            cur_seq.clear();
        } else {
            cur_seq += line;
        }
    }
    if (!cur_acc.empty()) seqs[cur_acc] = cur_seq;
    return seqs;
}

// Normalize a nucleotide sequence: uppercase, non-ACGT -> N.
// BLAST DB stores ambiguous bases as N, while NCBI efetch returns
// IUPAC ambiguity codes (R, Y, S, W, K, M, B, D, H, V).
// This normalization ensures a fair comparison.
static std::string normalize_seq(const std::string& seq) {
    std::string result = seq;
    for (auto& c : result) {
        c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') c = 'N';
    }
    return result;
}

// Build OutputHit entries for the target accessions (full sequence, forward strand).
static std::vector<OutputHit> make_hits() {
    std::vector<OutputHit> hits;
    for (int i = 0; i < NUM_TARGETS; i++) {
        if (!g_acc_info[i].found) continue;
        OutputHit h;
        h.query_id = "ssu_test";
        h.accession = TARGET_ACC[i];
        h.strand = '+';
        h.q_start = 0;
        h.q_end = g_acc_info[i].seq_length - 1;
        h.s_start = 0;
        h.s_end = g_acc_info[i].seq_length - 1;
        h.score = 100;
        h.volume = 0;
        hits.push_back(h);
    }
    return hits;
}

// ---- Setup ----

static bool setup_database() {
    g_db_dir = "/tmp/ikafssn_ssu_test";
    g_db_prefix = g_db_dir + "/SSU_eukaryote_rRNA";

    std::filesystem::create_directories(g_db_dir);

    // Check if already extracted (single- or multi-volume)
    bool have_db = std::filesystem::exists(g_db_prefix + ".nsq") ||
                   std::filesystem::exists(g_db_prefix + ".00.nsq");
    if (have_db) {
        std::fprintf(stderr, "SSU_eukaryote_rRNA already present, skipping download\n");
    } else {
        // Download
        std::string tar_path = g_db_dir + "/SSU_eukaryote_rRNA.tar.gz";
        std::string cmd = "curl -fSL --connect-timeout 30 --max-time 300 -o "
                          + tar_path + " " + SSU_URL;
        std::fprintf(stderr, "Downloading SSU_eukaryote_rRNA...\n");
        int ret = std::system(cmd.c_str());
        if (ret != 0) {
            std::fprintf(stderr, "Download failed (exit code %d)\n", ret);
            return false;
        }

        // Extract
        cmd = "tar xzf " + tar_path + " -C " + g_db_dir;
        std::fprintf(stderr, "Extracting...\n");
        ret = std::system(cmd.c_str());
        if (ret != 0) {
            std::fprintf(stderr, "Extraction failed (exit code %d)\n", ret);
            return false;
        }
    }

    // Open DB and find the target accessions
    auto vol_paths = BlastDbReader::find_volume_paths(g_db_prefix);
    if (vol_paths.empty()) {
        std::fprintf(stderr, "No volumes found for %s\n", g_db_prefix.c_str());
        return false;
    }
    std::fprintf(stderr, "Found %zu volume(s)\n", vol_paths.size());

    for (int i = 0; i < NUM_TARGETS; i++) g_acc_info[i].found = false;

    for (const auto& vpath : vol_paths) {
        BlastDbReader reader;
        if (!reader.open(vpath)) continue;

        uint32_t nseqs = reader.num_sequences();
        for (uint32_t oid = 0; oid < nseqs; oid++) {
            std::string acc = reader.get_accession(oid);
            for (int i = 0; i < NUM_TARGETS; i++) {
                if (acc == TARGET_ACC[i]) {
                    g_acc_info[i].seq_length = reader.seq_length(oid);
                    g_acc_info[i].found = true;
                    std::fprintf(stderr, "  Found %s (length=%u)\n",
                                 TARGET_ACC[i], g_acc_info[i].seq_length);
                }
            }
        }
    }

    return true;
}

// ---- Test 1: Local BLAST DB retrieval ----

static void test_local_retrieval() {
    std::fprintf(stderr, "-- test_local_retrieval\n");

    // All 3 accessions must exist in the DB
    for (int i = 0; i < NUM_TARGETS; i++) {
        if (!g_acc_info[i].found) {
            std::fprintf(stderr, "  %s NOT found in SSU_eukaryote_rRNA\n", TARGET_ACC[i]);
        }
        CHECK(g_acc_info[i].found);
        CHECK(g_acc_info[i].seq_length > 0);
    }

    auto hits = make_hits();
    CHECK_EQ(hits.size(), static_cast<size_t>(NUM_TARGETS));

    RetrieveOptions opts;
    opts.context = 0;

    std::ostringstream out;
    uint32_t retrieved = retrieve_local(hits, g_db_prefix, opts, out);
    CHECK_EQ(retrieved, static_cast<uint32_t>(NUM_TARGETS));

    auto seqs = parse_fasta_output(out.str());

    for (int i = 0; i < NUM_TARGETS; i++) {
        auto it = seqs.find(TARGET_ACC[i]);
        CHECK(it != seqs.end());
        if (it != seqs.end()) {
            // Sequence length must match the DB metadata
            CHECK_EQ(it->second.size(), static_cast<size_t>(g_acc_info[i].seq_length));
            // Must contain only valid bases (ACGTN)
            for (char c : it->second) {
                CHECK(c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N');
            }
            std::fprintf(stderr, "  %s: OK (%zu bp)\n",
                         TARGET_ACC[i], it->second.size());
        }
    }
}

// ---- Test 2: NCBI efetch remote retrieval ----

#ifdef IKAFSSN_ENABLE_REMOTE

// Saved for comparison in test 3
static std::string g_remote_fasta;

static void test_remote_retrieval() {
    std::fprintf(stderr, "-- test_remote_retrieval\n");

    auto hits = make_hits();
    CHECK_EQ(hits.size(), static_cast<size_t>(NUM_TARGETS));

    EfetchOptions opts;
    opts.context = 0;
    opts.batch_size = 100;
    opts.retries = 3;
    opts.timeout_sec = 60;
    opts.range_threshold = 100000;

    // Use API key from environment if available
    const char* env_key = std::getenv("NCBI_API_KEY");
    if (env_key) opts.api_key = env_key;

    std::ostringstream out;
    uint32_t retrieved = retrieve_remote(hits, opts, out);
    g_remote_fasta = out.str();

    CHECK_EQ(retrieved, static_cast<uint32_t>(NUM_TARGETS));

    auto seqs = parse_fasta_output(g_remote_fasta);

    for (int i = 0; i < NUM_TARGETS; i++) {
        auto it = seqs.find(TARGET_ACC[i]);
        CHECK(it != seqs.end());
        if (it != seqs.end()) {
            CHECK(it->second.size() > 0);
            std::fprintf(stderr, "  %s: OK (%zu bp)\n",
                         TARGET_ACC[i], it->second.size());
        }
    }
}

// ---- Test 3: Local vs remote comparison ----

static void test_local_vs_remote() {
    std::fprintf(stderr, "-- test_local_vs_remote\n");

    // Retrieve locally (fresh call)
    auto hits = make_hits();
    RetrieveOptions local_opts;
    std::ostringstream local_out;
    retrieve_local(hits, g_db_prefix, local_opts, local_out);

    auto local_seqs = parse_fasta_output(local_out.str());
    auto remote_seqs = parse_fasta_output(g_remote_fasta);

    for (int i = 0; i < NUM_TARGETS; i++) {
        auto local_it = local_seqs.find(TARGET_ACC[i]);
        auto remote_it = remote_seqs.find(TARGET_ACC[i]);

        if (local_it == local_seqs.end() || remote_it == remote_seqs.end()) {
            std::fprintf(stderr, "  %s: MISSING in local or remote output\n",
                         TARGET_ACC[i]);
            CHECK(false);
            continue;
        }

        // Normalize: BLAST DB maps ambiguous bases to N, efetch uses IUPAC codes
        std::string local_norm = normalize_seq(local_it->second);
        std::string remote_norm = normalize_seq(remote_it->second);

        CHECK_EQ(local_norm.size(), remote_norm.size());

        bool match = (local_norm == remote_norm);
        CHECK(match);

        if (match) {
            std::fprintf(stderr, "  %s: IDENTICAL (%zu bp)\n",
                         TARGET_ACC[i], local_norm.size());
        } else {
            // Report first mismatch position for debugging
            size_t min_len = std::min(local_norm.size(), remote_norm.size());
            for (size_t j = 0; j < min_len; j++) {
                if (local_norm[j] != remote_norm[j]) {
                    std::fprintf(stderr,
                        "  %s: MISMATCH at pos %zu: local='%c' remote='%c' "
                        "(local_raw='%c' remote_raw='%c')\n",
                        TARGET_ACC[i], j, local_norm[j], remote_norm[j],
                        local_it->second[j], remote_it->second[j]);
                    break;
                }
            }
            if (local_norm.size() != remote_norm.size()) {
                std::fprintf(stderr,
                    "  %s: LENGTH MISMATCH local=%zu remote=%zu\n",
                    TARGET_ACC[i], local_norm.size(), remote_norm.size());
            }
        }
    }
}

#endif // IKAFSSN_ENABLE_REMOTE

// ---- Main ----

int main() {
    if (!setup_database()) {
        std::fprintf(stderr,
            "Failed to set up SSU_eukaryote_rRNA database, skipping all tests\n");
        // Return 0 (skip) rather than 1 (fail) for network issues
        return 0;
    }

    test_local_retrieval();

#ifdef IKAFSSN_ENABLE_REMOTE
    test_remote_retrieval();
    test_local_vs_remote();
#else
    std::fprintf(stderr,
        "SKIPPED: remote and comparison tests (built without ENABLE_REMOTE_RETRIEVE)\n");
#endif

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
