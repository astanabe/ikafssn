#pragma once

// Shared test fixture for SSU_eukaryote_rRNA-based tests.
// Provides paths, constants, and helper functions.

#include "io/blastdb_reader.hpp"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <string>

#ifndef SOURCE_DIR
#define SOURCE_DIR "."
#endif

namespace ssu_fixture {

// ---- Paths ----

inline std::string ssu_db_prefix() {
    return std::string(SOURCE_DIR) + "/db/SSU_eukaryote_rRNA";
}

inline const char* DERIVED_DIR = "/tmp/ikafssn_ssu_test";

inline std::string ambig_db_prefix() {
    return std::string(DERIVED_DIR) + "/ssu_ambigdb";
}

inline std::string multivol_a_prefix() {
    return std::string(DERIVED_DIR) + "/ssu_multivol_a";
}

inline std::string multivol_b_prefix() {
    return std::string(DERIVED_DIR) + "/ssu_multivol_b";
}

inline std::string queries_path() {
    return std::string(DERIVED_DIR) + "/queries.fasta";
}

inline std::string seqidlist_path() {
    return std::string(DERIVED_DIR) + "/seqidlist.txt";
}

// ---- Target accession constants ----
// BlastDbReader::get_accession() returns unversioned accessions (NCBI convention).
// Use versioned form (ACC_FJ_V etc.) for blastdbcmd and seqidlist.

inline const char* ACC_FJ = "FJ876973";
inline const char* ACC_GQ = "GQ912721";
inline const char* ACC_DQ = "DQ235612";

// Versioned forms for external tools (blastdbcmd, seqidlist)
inline const char* ACC_FJ_V = "FJ876973.1";
inline const char* ACC_GQ_V = "GQ912721.1";
inline const char* ACC_DQ_V = "DQ235612.1";

// ---- Helpers ----

// Skip test (CTest interprets exit(0) as SKIP when no assertions fail).
[[noreturn]] inline void skip(const char* reason) {
    std::fprintf(stderr, "SKIP: %s\n", reason);
    std::exit(0);
}

// Verify the SSU DB is available; skip if not.
inline void check_ssu_available() {
    if (!std::filesystem::exists(ssu_db_prefix() + ".nsq") &&
        !std::filesystem::exists(ssu_db_prefix() + ".00.nsq")) {
        skip("SSU_eukaryote_rRNA not found at db/SSU_eukaryote_rRNA");
    }
}

// Verify derived test data is ready; skip if not.
inline void check_derived_data_ready() {
    if (!std::filesystem::exists(ambig_db_prefix() + ".nsq")) {
        skip("Derived test data not found. Run test/scripts/setup_ssu_testdata.sh first.");
    }
}

// Find OID by accession via linear scan. Returns UINT32_MAX if not found.
inline uint32_t find_oid_by_accession(ikafssn::BlastDbReader& db,
                                      const std::string& acc) {
    uint32_t n = db.num_sequences();
    for (uint32_t oid = 0; oid < n; oid++) {
        std::string a = db.get_accession(oid);
        if (a == acc) return oid;
    }
    return UINT32_MAX;
}

// Extract a subsequence string from the DB for a given accession and range.
// Returns empty string on failure.
inline std::string extract_subsequence(ikafssn::BlastDbReader& db,
                                       const std::string& acc,
                                       uint32_t start, uint32_t end) {
    uint32_t oid = find_oid_by_accession(db, acc);
    if (oid == UINT32_MAX) return {};
    std::string full_seq = db.get_sequence(oid);
    if (end > full_seq.size()) end = static_cast<uint32_t>(full_seq.size());
    if (start >= end) return {};
    return full_seq.substr(start, end - start);
}

} // namespace ssu_fixture
