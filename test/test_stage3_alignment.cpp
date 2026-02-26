#include "test_util.hpp"

#include "search/stage3_alignment.hpp"
#include "io/fasta_reader.hpp"
#include "io/result_writer.hpp"
#include "util/logger.hpp"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include <parasail.h>

using namespace ikafssn;

// ---- Parasail direct tests (no BLAST DB needed) ----

static void test_parasail_exact_match() {
    std::fprintf(stderr, "-- test_parasail_exact_match\n");

    // Align identical sequences using parasail directly
    const char* seq = "ACGTACGTACGTACGTACGT";
    int len = static_cast<int>(std::strlen(seq));
    const parasail_matrix_t* matrix = parasail_matrix_lookup("nuc44");

    parasail_result_t* result = parasail_sg_trace_striped_sat(
        seq, len, seq, len, 10, 1, matrix);

    CHECK(result != nullptr);
    CHECK(result->score > 0);

    // Get CIGAR
    parasail_cigar_t* cigar = parasail_result_get_cigar(
        result, seq, len, seq, len, matrix);
    CHECK(cigar != nullptr);
    CHECK(cigar->len > 0);

    // Decode CIGAR - should be all '=' for exact match
    char* cigar_str = parasail_cigar_decode(cigar);
    CHECK(cigar_str != nullptr);
    std::string cs(cigar_str);
    // Should be "20=" for exact match
    CHECK(cs == "20=");

    free(cigar_str);
    parasail_cigar_free(cigar);
    parasail_result_free(result);
}

static void test_parasail_with_mismatches() {
    std::fprintf(stderr, "-- test_parasail_with_mismatches\n");

    const char* query = "ACGTACGTACGTACGTACGT";
    const char* subj  = "ACGTACTTACGTACGTACGT";
    //                          ^  mismatch at pos 6 (G->T)
    int qlen = static_cast<int>(std::strlen(query));
    int slen = static_cast<int>(std::strlen(subj));
    const parasail_matrix_t* matrix = parasail_matrix_lookup("nuc44");

    parasail_result_t* result = parasail_sg_trace_striped_sat(
        query, qlen, subj, slen, 10, 1, matrix);
    CHECK(result != nullptr);

    parasail_cigar_t* cigar = parasail_result_get_cigar(
        result, query, qlen, subj, slen, matrix);
    CHECK(cigar != nullptr);

    // Walk CIGAR to count matches and mismatches
    uint32_t nident = 0, nmismatch = 0;
    for (int i = 0; i < cigar->len; i++) {
        char op = parasail_cigar_decode_op(cigar->seq[i]);
        uint32_t len = parasail_cigar_decode_len(cigar->seq[i]);
        if (op == '=') nident += len;
        else if (op == 'X') nmismatch += len;
    }
    CHECK(nident == 19);
    CHECK(nmismatch == 1);

    parasail_cigar_free(cigar);
    parasail_result_free(result);
}

static void test_parasail_with_gaps() {
    std::fprintf(stderr, "-- test_parasail_with_gaps\n");

    // Query has a deletion relative to subject
    const char* query = "ACGTACGTACGT";       // 12bp
    const char* subj  = "ACGTAAACGTACGT";     // 14bp (extra "AA" inserted)
    int qlen = static_cast<int>(std::strlen(query));
    int slen = static_cast<int>(std::strlen(subj));
    const parasail_matrix_t* matrix = parasail_matrix_lookup("nuc44");

    parasail_result_t* result = parasail_sg_trace_striped_sat(
        query, qlen, subj, slen, 10, 1, matrix);
    CHECK(result != nullptr);

    parasail_cigar_t* cigar = parasail_result_get_cigar(
        result, query, qlen, subj, slen, matrix);
    CHECK(cigar != nullptr);

    // Decode and check that CIGAR contains I or D
    char* cigar_str = parasail_cigar_decode(cigar);
    std::string cs(cigar_str);
    // Should have D (deletion in query = insertion in ref) or I
    CHECK(cs.find('D') != std::string::npos || cs.find('I') != std::string::npos);

    free(cigar_str);
    parasail_cigar_free(cigar);
    parasail_result_free(result);
}

static void test_parasail_score_only() {
    std::fprintf(stderr, "-- test_parasail_score_only\n");

    const char* seq = "ACGTACGTACGTACGTACGT";
    int len = static_cast<int>(std::strlen(seq));
    const parasail_matrix_t* matrix = parasail_matrix_lookup("nuc44");

    // Score-only (no trace)
    parasail_result_t* result = parasail_sg_striped_sat(
        seq, len, seq, len, 10, 1, matrix);
    CHECK(result != nullptr);
    CHECK(result->score > 0);
    CHECK(result->end_query >= 0);
    CHECK(result->end_ref >= 0);

    parasail_result_free(result);
}

static void test_parasail_profile() {
    std::fprintf(stderr, "-- test_parasail_profile\n");

    const char* query = "ACGTACGTACGTACGTACGT";
    const char* subj  = "ACGTACGTACGTACGTACGT";
    int qlen = static_cast<int>(std::strlen(query));
    int slen = static_cast<int>(std::strlen(subj));
    const parasail_matrix_t* matrix = parasail_matrix_lookup("nuc44");

    parasail_profile_t* profile = parasail_profile_create_sat(query, qlen, matrix);
    CHECK(profile != nullptr);

    // Score-only with profile
    parasail_result_t* result = parasail_sg_striped_profile_sat(
        profile, subj, slen, 10, 1);
    CHECK(result != nullptr);
    CHECK(result->score > 0);

    parasail_result_free(result);

    // Trace with profile
    parasail_result_t* result2 = parasail_sg_trace_striped_profile_sat(
        profile, subj, slen, 10, 1);
    CHECK(result2 != nullptr);
    CHECK(result2->score > 0);

    parasail_result_free(result2);
    parasail_profile_free(profile);
}

static void test_reverse_complement_alignment() {
    std::fprintf(stderr, "-- test_reverse_complement_alignment\n");

    // Test that RC of a sequence aligns perfectly to itself
    const char* fwd = "ACGTACGTACGT";
    int fwd_len = static_cast<int>(std::strlen(fwd));
    const parasail_matrix_t* matrix = parasail_matrix_lookup("nuc44");

    // Manual RC
    std::string rc(fwd);
    std::reverse(rc.begin(), rc.end());
    for (auto& c : rc) {
        switch (c) {
            case 'A': c = 'T'; break;
            case 'T': c = 'A'; break;
            case 'C': c = 'G'; break;
            case 'G': c = 'C'; break;
        }
    }

    // Align RC to itself should be perfect
    parasail_result_t* result = parasail_sg_trace_striped_sat(
        rc.c_str(), static_cast<int>(rc.size()),
        rc.c_str(), static_cast<int>(rc.size()),
        10, 1, matrix);
    CHECK(result != nullptr);
    CHECK(result->score > 0);

    parasail_cigar_t* cigar = parasail_result_get_cigar(
        result, rc.c_str(), static_cast<int>(rc.size()),
        rc.c_str(), static_cast<int>(rc.size()), matrix);
    char* cigar_str = parasail_cigar_decode(cigar);
    std::string cs(cigar_str);
    CHECK(cs == "12="); // perfect match

    free(cigar_str);
    parasail_cigar_free(cigar);
    parasail_result_free(result);
}

static void test_pident_calculation() {
    std::fprintf(stderr, "-- test_pident_calculation\n");

    // 18 matches + 2 mismatches = 20 bases, pident = 90%
    const char* query = "ACGTACGTACGTACGTACGT";
    const char* subj  = "ACTTACGTACGTACGTACTT";
    //                     ^                ^^  mismatches
    int qlen = static_cast<int>(std::strlen(query));
    int slen = static_cast<int>(std::strlen(subj));
    const parasail_matrix_t* matrix = parasail_matrix_lookup("nuc44");

    parasail_result_t* result = parasail_sg_trace_striped_sat(
        query, qlen, subj, slen, 10, 1, matrix);
    CHECK(result != nullptr);

    parasail_cigar_t* cigar = parasail_result_get_cigar(
        result, query, qlen, subj, slen, matrix);
    CHECK(cigar != nullptr);

    uint32_t nident = 0, nmismatch = 0, aln_len = 0;
    for (int i = 0; i < cigar->len; i++) {
        char op = parasail_cigar_decode_op(cigar->seq[i]);
        uint32_t len = parasail_cigar_decode_len(cigar->seq[i]);
        aln_len += len;
        if (op == '=') nident += len;
        else if (op == 'X') nmismatch += len;
    }

    double pident = (aln_len > 0) ? 100.0 * nident / aln_len : 0.0;
    CHECK(nident + nmismatch == 20);
    CHECK(pident > 85.0 && pident <= 100.0);

    parasail_cigar_free(cigar);
    parasail_result_free(result);
}

int main() {
    test_parasail_exact_match();
    test_parasail_with_mismatches();
    test_parasail_with_gaps();
    test_parasail_score_only();
    test_parasail_profile();
    test_reverse_complement_alignment();
    test_pident_calculation();

    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
