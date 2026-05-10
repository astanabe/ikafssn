// unit test: src/search/result_dedup.cpp.
//
// The Stage 2 dedup collapses orchestrator hits whose parent-relative
// (sstart, send) — derived by adding fragment_start - 1 — coincide,
// keeping a single canonical row per (query, volume, parent, strand,
// parent_sstart, parent_send, chainscore) tuple.  The Stage 3 dedup is
// the post-alignment counterpart, keyed on
// (qseqid, sseqid, sstrand, send, alnscore) and skipping any
// skip_reason != 0 rows.
//
// We feed each helper synthetic input that mimics the duplicate that two
// adjacent fragments would produce (different fragment-relative s_start /
// s_end, but the same parent-relative coordinates after the
// fragment_start shift) and assert exactly one row survives.

#include "test_util.hpp"

#include "index/ksx_writer.hpp"
#include "index/ksx_reader.hpp"
#include "io/result_writer.hpp"
#include "search/parallel_search.hpp"
#include "search/result_dedup.hpp"

#include <cstdio>
#include <string>
#include <vector>

using namespace ikafssn;

namespace {

OrchestratorHit make_hit(size_t query_idx, size_t vol_idx, uint16_t volume_index,
                        uint32_t seq_id, bool is_reverse,
                        uint32_t s_start, uint32_t s_end,
                        uint32_t chainscore) {
    OrchestratorHit h{};
    h.query_idx    = query_idx;
    h.volume_idx   = vol_idx;
    h.volume_index = volume_index;
    h.cr.seq_id     = seq_id;
    h.cr.is_reverse = is_reverse;
    h.cr.s_start    = s_start;
    h.cr.s_end      = s_end;
    h.cr.q_start    = 1;
    h.cr.q_end      = 100;
    h.cr.chainscore = chainscore;
    return h;
}

OutputHit make_output_hit(const std::string& q, const std::string& s,
                          char strand, uint32_t sstart, uint32_t send,
                          int32_t alnscore, uint8_t skip_reason = 0) {
    OutputHit h;
    h.qseqid  = q;
    h.sseqid  = s;
    h.sstrand = strand;
    h.qstart  = 1;
    h.qend    = 100;
    h.sstart  = sstart;
    h.send    = send;
    h.alnscore    = alnscore;
    h.skip_reason = skip_reason;
    return h;
}

// Build a small two-fragment .ksx so that Stage 2 dedup has a real
// KsxReader to query for parent_index() / fragment_start().  Two
// fragments of one parent with the LONGCHR_1-style overlap.
std::string write_test_ksx() {
    static const char* path = "/tmp/test_result_dedup.ksx";
    KsxWriter writer;
    writer.set_min_seq_length(64);
    writer.set_min_length_split(1500);
    writer.set_overlap_length(200);
    uint32_t p0 = writer.add_parent(0, 5171, "LONGCHR_1");
    writer.add_fragment(p0, 1,    1791);  // SeqId 0
    writer.add_fragment(p0, 1592, 3381);  // SeqId 1
    writer.add_fragment(p0, 3182, 4771);  // SeqId 2
    uint32_t p1 = writer.add_parent(1, 800, "SHORT");
    writer.add_fragment(p1, 1, 800);      // SeqId 3
    if (!writer.write(path)) {
        std::fprintf(stderr, "FATAL: failed to write %s\n", path);
        std::exit(1);
    }
    return path;
}

}  // namespace

// 1. Two adjacent fragments of the same parent emit chains whose
// fragment-relative coordinates land on the same parent-relative range.
// Stage 2 dedup collapses them to a single hit.
static void test_stage2_collapses_adjacent_fragment_dups() {
    std::fprintf(stderr, "-- test_stage2_collapses_adjacent_fragment_dups\n");

    std::string ksx_path = write_test_ksx();
    KsxReader ksx;
    CHECK(ksx.open(ksx_path));

    // Hit A: produced from SeqId 0 (fragment_start=1) at fragment-relative
    // s_start=1700, s_end=1789  -> parent-relative [1700, 1789].
    // Hit B: produced from SeqId 1 (fragment_start=1592) at fragment-
    // relative s_start=109, s_end=198  -> parent-relative [1700, 1789].
    // Same query, same volume, same strand, same chainscore: dedup must
    // collapse these to a single row.
    std::vector<OrchestratorHit> hits = {
        make_hit(/*qi*/ 0, /*vi*/ 0, /*vol*/ 0, /*sid*/ 0,
                 /*rev*/ false, /*ss*/ 1700, /*se*/ 1789, /*cs*/ 42),
        make_hit(0, 0, 0, 1, false, 109, 198, 42),
    };

    std::vector<const KsxReader*> ksx_per_volume = { &ksx };
    dedup_stage2_orchestrator_hits(hits, ksx_per_volume);
    CHECK_EQ(hits.size(), 1u);

    // Negative control: a hit on the OTHER parent (SeqId 3) is NOT a
    // duplicate of the LONGCHR_1 hit and must survive alongside it.
    hits = {
        make_hit(0, 0, 0, 0, false, 1700, 1789, 42),
        make_hit(0, 0, 0, 1, false, 109,  198,  42),
        make_hit(0, 0, 0, 3, false, 100,  189,  42),
    };
    dedup_stage2_orchestrator_hits(hits, ksx_per_volume);
    CHECK_EQ(hits.size(), 2u);

    // Different chainscore -> different key, both rows kept (the dedup is
    // intentionally score-aware so the per-fragment best chain stays).
    hits = {
        make_hit(0, 0, 0, 0, false, 1700, 1789, 40),
        make_hit(0, 0, 0, 1, false, 109,  198,  42),
    };
    dedup_stage2_orchestrator_hits(hits, ksx_per_volume);
    CHECK_EQ(hits.size(), 2u);

    // Forward and reverse strands of the same parent-relative range are
    // separate keys -> both kept.
    hits = {
        make_hit(0, 0, 0, 0, /*rev*/ false, 1700, 1789, 42),
        make_hit(0, 0, 0, 0, /*rev*/ true,  1700, 1789, 42),
    };
    dedup_stage2_orchestrator_hits(hits, ksx_per_volume);
    CHECK_EQ(hits.size(), 2u);

    ksx.close();
    std::remove(ksx_path.c_str());
}

// 2. Stage 3 dedup keys: (qseqid, sseqid, sstrand, send, alnscore).
// Two rows with the same key but different qstart / sstart / coverscore
// (the kind of variance two adjacent fragments of one parent would
// produce post-Stage-3) collapse to one.
static void test_stage3_collapses_post_align_dups() {
    std::fprintf(stderr, "-- test_stage3_collapses_post_align_dups\n");

    std::vector<OutputHit> hits = {
        make_output_hit("q1", "LONGCHR_1", '+', /*ss*/ 1700, /*se*/ 1789, /*aln*/ 178),
        make_output_hit("q1", "LONGCHR_1", '+', /*ss*/ 1701, /*se*/ 1789, /*aln*/ 178),
    };
    // Mutate qstart on one row so it's clear those fields don't drive the key.
    hits[1].qstart = 2;
    dedup_stage3_output_hits(hits);
    CHECK_EQ(hits.size(), 1u);

    // Different alnscore -> different key, both rows survive.
    hits = {
        make_output_hit("q1", "LONGCHR_1", '+', 1700, 1789, 178),
        make_output_hit("q1", "LONGCHR_1", '+', 1700, 1789, 175),
    };
    dedup_stage3_output_hits(hits);
    CHECK_EQ(hits.size(), 2u);

    // Different sseqid (parent accession) -> not a duplicate.
    hits = {
        make_output_hit("q1", "LONGCHR_1", '+', 1700, 1789, 178),
        make_output_hit("q1", "OTHER",     '+', 1700, 1789, 178),
    };
    dedup_stage3_output_hits(hits);
    CHECK_EQ(hits.size(), 2u);
}

// 3. Skip-marker rows (skip_reason != 0) are preserved verbatim and not
// run through the dedup, so two identical skip rows stay identical
// (callers may emit one skip per query and rely on that).
static void test_stage3_preserves_skip_rows() {
    std::fprintf(stderr, "-- test_stage3_preserves_skip_rows\n");
    std::vector<OutputHit> hits = {
        make_output_hit("q1", "LONGCHR_1", '+', 1700, 1789, 178),
        make_output_hit("q1", "LONGCHR_1", '+', 1700, 1789, 178),
        make_output_hit("q2", "",          '+', 0, 0, 0, /*skip*/ 5),
        make_output_hit("q2", "",          '+', 0, 0, 0, /*skip*/ 5),
    };
    dedup_stage3_output_hits(hits);
    // Real hits collapse to 1; skip rows are passed through unchanged (2).
    CHECK_EQ(hits.size(), 3u);

    // Real hits land before skip markers (Stage 3 dedup splits the vector
    // and re-appends skip rows at the end).
    size_t real_count = 0, skip_count = 0;
    for (auto& h : hits) {
        if (h.skip_reason == 0) real_count++;
        else                    skip_count++;
    }
    CHECK_EQ(real_count, 1u);
    CHECK_EQ(skip_count, 2u);
}

// 4. Stage 2 dedup with no input / single input is a no-op and must not
// touch the supplied ksx_per_volume vector.
static void test_stage2_trivial_inputs() {
    std::fprintf(stderr, "-- test_stage2_trivial_inputs\n");
    std::vector<const KsxReader*> ksx_per_volume = { nullptr };

    std::vector<OrchestratorHit> empty;
    dedup_stage2_orchestrator_hits(empty, ksx_per_volume);
    CHECK_EQ(empty.size(), 0u);

    std::vector<OrchestratorHit> singleton = {
        make_hit(0, 0, 0, 0, false, 1, 100, 42),
    };
    dedup_stage2_orchestrator_hits(singleton, ksx_per_volume);
    CHECK_EQ(singleton.size(), 1u);
}

int main() {
    test_stage2_collapses_adjacent_fragment_dups();
    test_stage3_collapses_post_align_dups();
    test_stage3_preserves_skip_rows();
    test_stage2_trivial_inputs();
    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
