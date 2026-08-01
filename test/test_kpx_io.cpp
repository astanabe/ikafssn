#include "test_util.hpp"
#include "kpx_writer.hpp"
#include "index/kpx_reader.hpp"
#include "core/config.hpp"
#include "index/pfd_codec.hpp"

#include <algorithm>
#include <cstdio>
#include <vector>
#include <string>

using namespace ikafssn;

static const char* TEST_FILE = "/tmp/test_ikafssn.kpx";

// Decode a .kpx posting list via the candidate-set API.  The writer adds
// postings only for `candidates`, so the sorted unique candidate list is
// also the k-mer's .kix decoded distinct seq_id array.  Returns the
// positions concatenated in candidate order.
static std::vector<uint32_t> decode_pos_postings(
    const uint8_t* data, const uint8_t* end,
    const std::vector<uint32_t>& candidates) {
    std::vector<uint32_t> distinct_sorted = candidates;
    std::sort(distinct_sorted.begin(), distinct_sorted.end());
    distinct_sorted.erase(std::unique(distinct_sorted.begin(),
                                       distinct_sorted.end()),
                          distinct_sorted.end());

    pfd::PosDecodeScratch scratch;
    if (!pfd::open_stream_kpx_for_candidates(
            data, static_cast<size_t>(end - data),
            distinct_sorted.data(), distinct_sorted.size(),
            distinct_sorted.data(), distinct_sorted.size(),
            scratch)) {
        return {};
    }
    return scratch.out_positions;
}

static void test_single_seq() {
    int k = 7;
    uint32_t ts = table_size(k);

    // Single sequence, multiple positions for one k-mer (occ>=2 path).
    std::vector<KpxWriter::PostingEntry> entries = {
        {5, 10}, {5, 20}, {5, 30}, {5, 100}
    };
    std::vector<uint32_t> candidates = {5};

    {
        KpxWriter writer(k);
        for (uint32_t i = 0; i < ts; i++) {
            if (i == 42) {
                writer.add_posting_list(i, entries);
            } else {
                writer.add_posting_list(i, {});
            }
        }
        CHECK(writer.write(TEST_FILE));
    }

    {
        KpxReader reader;
        CHECK(reader.open(TEST_FILE));
        CHECK_EQ(reader.k(), k);
        CHECK_EQ(reader.total_position_count(), 4u);

        auto positions = decode_pos_postings(
            reader.posting_file() + reader.pos_offset(42),
            reader.posting_file() + reader.posting_file_size(),
            candidates);

        CHECK_EQ(positions.size(), 4u);
        CHECK_EQ(positions[0], 10u);
        CHECK_EQ(positions[1], 20u);
        CHECK_EQ(positions[2], 30u);
        CHECK_EQ(positions[3], 100u);

        reader.close();
    }

    std::remove(TEST_FILE);
}

// All occ=1 clusters → entire posting list goes through the short_occ1
// sub-bucket (no u8 occ_count[] bytes, no partition groups).
static void test_all_occ1() {
    int k = 5;
    uint32_t ts = table_size(k);

    std::vector<KpxWriter::PostingEntry> entries;
    std::vector<uint32_t> candidates;
    for (uint32_t sid = 0; sid < 50; sid++) {
        entries.push_back({sid, 100u + sid});
        candidates.push_back(sid);
    }

    {
        KpxWriter writer(k);
        for (uint32_t i = 0; i < ts; i++) {
            if (i == 17) writer.add_posting_list(i, entries);
            else         writer.add_posting_list(i, {});
        }
        CHECK(writer.write(TEST_FILE));
    }

    {
        KpxReader reader;
        CHECK(reader.open(TEST_FILE));
        CHECK_EQ(reader.total_position_count(), 50u);

        auto positions = decode_pos_postings(
            reader.posting_file() + reader.pos_offset(17),
            reader.posting_file() + reader.posting_file_size(),
            candidates);
        CHECK_EQ(positions.size(), 50u);
        for (uint32_t i = 0; i < 50; i++) {
            CHECK_EQ(positions[i], 100u + i);
        }
        reader.close();
    }
    std::remove(TEST_FILE);
}

// Mixed occ=1 / occ>=2 (no partition since default threshold = 8 and all
// occ counts <= 8).  Exercises both short sub-buckets.
static void test_mixed_short_buckets() {
    int k = 5;
    uint32_t ts = table_size(k);

    std::vector<KpxWriter::PostingEntry> entries = {
        {1, 5},                 // occ=1
        {3, 10}, {3, 20},       // occ=2 → short_occ_ge2
        {7, 99},                // occ=1
        {9, 1}, {9, 2}, {9, 3}, // occ=3
        {11, 50},               // occ=1
    };
    std::vector<uint32_t> candidates = {1, 3, 7, 9, 11};

    {
        KpxWriter writer(k);
        for (uint32_t i = 0; i < ts; i++) {
            if (i == 5) writer.add_posting_list(i, entries);
            else        writer.add_posting_list(i, {});
        }
        CHECK(writer.write(TEST_FILE));
    }

    {
        KpxReader reader;
        CHECK(reader.open(TEST_FILE));
        CHECK_EQ(reader.total_position_count(), 8u);

        auto positions = decode_pos_postings(
            reader.posting_file() + reader.pos_offset(5),
            reader.posting_file() + reader.posting_file_size(),
            candidates);

        CHECK_EQ(positions.size(), 8u);
        CHECK_EQ(positions[0], 5u);    // sid=1 occ=1
        CHECK_EQ(positions[1], 10u);   // sid=3
        CHECK_EQ(positions[2], 20u);
        CHECK_EQ(positions[3], 99u);   // sid=7
        CHECK_EQ(positions[4], 1u);    // sid=9
        CHECK_EQ(positions[5], 2u);
        CHECK_EQ(positions[6], 3u);
        CHECK_EQ(positions[7], 50u);   // sid=11
        reader.close();
    }
    std::remove(TEST_FILE);
}

// Partition-only case: every cluster exceeds the threshold.  Hits the
// partition path with 0 short_occ1 / short_occ_ge2 sub-buckets.
static void test_partition_only() {
    int k = 5;
    uint32_t ts = table_size(k);

    std::vector<KpxWriter::PostingEntry> entries;
    std::vector<uint32_t> candidates;
    for (uint32_t sid : {3u, 4u, 5u}) {
        for (uint32_t pos = 0; pos < 12; pos++) {
            entries.push_back({sid, 1000u + sid * 100 + pos});
        }
        candidates.push_back(sid);
    }

    {
        KpxWriter writer(k, /*freq_threshold_part=*/8);  // 12 > 8
        for (uint32_t i = 0; i < ts; i++) {
            if (i == 21) writer.add_posting_list(i, entries);
            else         writer.add_posting_list(i, {});
        }
        CHECK(writer.write(TEST_FILE));
    }

    {
        KpxReader reader;
        CHECK(reader.open(TEST_FILE));
        CHECK_EQ(reader.total_position_count(), 36u);

        auto positions = decode_pos_postings(
            reader.posting_file() + reader.pos_offset(21),
            reader.posting_file() + reader.posting_file_size(),
            candidates);

        CHECK_EQ(positions.size(), 36u);
        for (size_t j = 0; j < entries.size(); j++) {
            CHECK_EQ(positions[j], entries[j].pos);
        }
        reader.close();
    }
    std::remove(TEST_FILE);
}

// FOR-block stream tail size sweep: 0, 1, 127 elements in the tail
// (alongside one full block of 128).  Single partition group with
// occurrence counts up to 255 (the u8 cap inherited from the dedup
// kernel).  For larger sweeps the short_occ1 sub-bucket is used in
// test_short1_full_blocks below.
static void test_partition_tail_sizes() {
    int k = 5;
    uint32_t ts = table_size(k);

    auto build_and_check = [&](uint32_t total_count) {
        std::vector<KpxWriter::PostingEntry> entries;
        entries.reserve(total_count);
        for (uint32_t pos = 0; pos < total_count; pos++) {
            entries.push_back({42u, 5000u + pos * 7u});
        }
        std::vector<uint32_t> candidates = {42};

        KpxWriter writer(k, /*freq_threshold_part=*/8);  // total_count > 8
        for (uint32_t i = 0; i < ts; i++) {
            if (i == 9) writer.add_posting_list(i, entries);
            else        writer.add_posting_list(i, {});
        }
        CHECK(writer.write(TEST_FILE));

        KpxReader reader;
        CHECK(reader.open(TEST_FILE));
        CHECK_EQ(reader.total_position_count(), total_count);

        auto positions = decode_pos_postings(
            reader.posting_file() + reader.pos_offset(9),
            reader.posting_file() + reader.posting_file_size(),
            candidates);
        CHECK_EQ(positions.size(), static_cast<size_t>(total_count));
        for (uint32_t i = 0; i < total_count; i++) {
            CHECK_EQ(positions[i], 5000u + i * 7u);
        }
        reader.close();
        std::remove(TEST_FILE);
    };

    build_and_check(128);   // 1 full block, tail_count = 0
    build_and_check(129);   // 1 full block, tail_count = 1
    build_and_check(255);   // 1 full block, tail_count = 127
}

// Regression: a single (k-mer, seq_id) cluster with > 255 occurrences.
// occ_count must stay wide through the dedup → encoder pipeline: a u8
// would saturate the run at 255, drop the rest of the positions, and
// corrupt the encoder's pos_cursor walk.
static void test_partition_above_255() {
    int k = 5;
    uint32_t ts = table_size(k);

    const uint32_t occ = 1000;  // well above the short bucket's u8 cap
    std::vector<KpxWriter::PostingEntry> entries;
    entries.reserve(occ);
    for (uint32_t pos = 0; pos < occ; pos++) {
        entries.push_back({77u, 100u + pos * 11u});
    }
    std::vector<uint32_t> candidates = {77};

    {
        KpxWriter writer(k, /*freq_threshold_part=*/8);  // occ=1000 → partition
        for (uint32_t i = 0; i < ts; i++) {
            if (i == 19) writer.add_posting_list(i, entries);
            else         writer.add_posting_list(i, {});
        }
        CHECK(writer.write(TEST_FILE));
    }

    {
        KpxReader reader;
        CHECK(reader.open(TEST_FILE));
        CHECK_EQ(reader.total_position_count(), occ);

        auto positions = decode_pos_postings(
            reader.posting_file() + reader.pos_offset(19),
            reader.posting_file() + reader.posting_file_size(),
            candidates);
        CHECK_EQ(positions.size(), static_cast<size_t>(occ));
        for (uint32_t i = 0; i < occ; i++) {
            CHECK_EQ(positions[i], 100u + i * 11u);
        }
        reader.close();
    }
    std::remove(TEST_FILE);
}

// 2+ full FOR blocks via the short_occ1 sub-bucket: each distinct sid
// contributes exactly 1 position, so the per-(kmer, seq_id) u8 occ cap
// does not bound the FOR stream length.
static void test_short1_full_blocks() {
    int k = 5;
    uint32_t ts = table_size(k);

    const uint32_t n = 256;  // 2 full FOR blocks, tail_count = 0
    std::vector<KpxWriter::PostingEntry> entries;
    std::vector<uint32_t> candidates;
    entries.reserve(n);
    candidates.reserve(n);
    for (uint32_t sid = 0; sid < n; sid++) {
        entries.push_back({sid, 100u + sid * 3u});
        candidates.push_back(sid);
    }

    {
        KpxWriter writer(k);
        for (uint32_t i = 0; i < ts; i++) {
            if (i == 11) writer.add_posting_list(i, entries);
            else         writer.add_posting_list(i, {});
        }
        CHECK(writer.write(TEST_FILE));
    }

    {
        KpxReader reader;
        CHECK(reader.open(TEST_FILE));
        CHECK_EQ(reader.total_position_count(), n);

        auto positions = decode_pos_postings(
            reader.posting_file() + reader.pos_offset(11),
            reader.posting_file() + reader.posting_file_size(),
            candidates);
        CHECK_EQ(positions.size(), static_cast<size_t>(n));
        for (uint32_t i = 0; i < n; i++) {
            CHECK_EQ(positions[i], 100u + i * 3u);
        }
        reader.close();
    }
    std::remove(TEST_FILE);
}

// distinct_count == 0 empty posting list (no k-mer matched, but writer still
// emits the all-zero blob so offsets stay valid).
static void test_empty_posting() {
    int k = 5;
    uint32_t ts = table_size(k);

    {
        KpxWriter writer(k);
        for (uint32_t i = 0; i < ts; i++) {
            writer.add_posting_list(i, {});
        }
        CHECK(writer.write(TEST_FILE));
    }

    {
        KpxReader reader;
        CHECK(reader.open(TEST_FILE));
        CHECK_EQ(reader.total_position_count(), 0u);

        auto positions = decode_pos_postings(
            reader.posting_file() + reader.pos_offset(0),
            reader.posting_file() + reader.posting_file_size(),
            {0});
        CHECK_EQ(positions.size(), 0u);
        reader.close();
    }
    std::remove(TEST_FILE);
}

// All candidates miss: posting list has data, but caller-supplied candidate
// set does not include any of its sids → every per-candidate bucket
// must come back empty.
static void test_all_candidates_miss() {
    int k = 5;
    uint32_t ts = table_size(k);

    std::vector<KpxWriter::PostingEntry> entries = {
        {1, 10}, {2, 20}, {3, 30}, {3, 31},
    };

    {
        KpxWriter writer(k);
        for (uint32_t i = 0; i < ts; i++) {
            if (i == 99) writer.add_posting_list(i, entries);
            else         writer.add_posting_list(i, {});
        }
        CHECK(writer.write(TEST_FILE));
    }

    {
        KpxReader reader;
        CHECK(reader.open(TEST_FILE));

        // The actual distinct sids in the posting list are {1, 2, 3}; the
        // .kix decoded array must reflect that.  Candidates {7, 8, 9}
        // intentionally have no overlap.
        std::vector<uint32_t> kix_decoded = {1, 2, 3};
        std::vector<uint32_t> candidates  = {7, 8, 9};
        pfd::PosDecodeScratch scratch;
        const uint64_t off = reader.pos_offset(99);
        CHECK(pfd::open_stream_kpx_for_candidates(
            reader.posting_file() + off,
            static_cast<size_t>(reader.posting_file_size() - off),
            kix_decoded.data(), kix_decoded.size(),
            candidates.data(), candidates.size(),
            scratch));
        CHECK_EQ(scratch.out_candidate_idx.size(), 0u);
        CHECK_EQ(scratch.out_positions.size(), 0u);
        reader.close();
    }
    std::remove(TEST_FILE);
}

// Candidate subset decode: only some of the posting list's distinct sids are
// candidates, and only some of the candidates are in the posting list.  The
// CSR must list the matched candidates by their index in the candidate array,
// which is what Stage 2A uses to address its per-candidate hit buckets.
static void test_candidate_subset_csr() {
    int k = 5;
    uint32_t ts = table_size(k);

    // sid 1 -> occ=1, sid 3 -> occ>=2, sid 5 -> occ=1, sid 7 -> partition.
    std::vector<KpxWriter::PostingEntry> entries = {
        {1, 10},
        {3, 20}, {3, 21}, {3, 22},
        {5, 30},
    };
    for (uint32_t j = 0; j < 12; j++) entries.push_back({7, 40 + j});  // 12 > 8

    {
        KpxWriter writer(k, /*freq_threshold_part=*/8);
        for (uint32_t i = 0; i < ts; i++) {
            if (i == 77) writer.add_posting_list(i, entries);
            else         writer.add_posting_list(i, {});
        }
        CHECK(writer.write(TEST_FILE));
    }

    {
        KpxReader reader;
        CHECK(reader.open(TEST_FILE));

        const std::vector<uint32_t> kix_decoded = {1, 3, 5, 7};
        // sid 1 is in the posting list but not a candidate; sids 0, 4, 9 are
        // candidates that the posting list does not contain.
        const std::vector<uint32_t> candidates = {0, 3, 4, 5, 7, 9};

        pfd::PosDecodeScratch scratch;
        const uint64_t off = reader.pos_offset(77);
        CHECK(pfd::open_stream_kpx_for_candidates(
            reader.posting_file() + off,
            static_cast<size_t>(reader.posting_file_size() - off),
            kix_decoded.data(), kix_decoded.size(),
            candidates.data(), candidates.size(),
            scratch));

        const std::vector<uint32_t> want_idx = {1, 3, 4};
        const std::vector<uint32_t> want_off = {0, 3, 4, 16};
        CHECK(scratch.out_candidate_idx == want_idx);
        CHECK(scratch.out_offsets == want_off);
        CHECK_EQ(scratch.out_positions.size(), 16u);

        // Positions of each matched candidate, addressed through the CSR.
        const std::vector<std::vector<uint32_t>> want_pos = {
            {20, 21, 22},
            {30},
            {40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51},
        };
        for (size_t m = 0; m < want_idx.size(); m++) {
            std::vector<uint32_t> got(
                scratch.out_positions.begin() + scratch.out_offsets[m],
                scratch.out_positions.begin() + scratch.out_offsets[m + 1]);
            CHECK(got == want_pos[m]);
        }
        reader.close();
    }
    std::remove(TEST_FILE);
}

// A posting list no candidate hits must be abandoned right after the
// selection pass: every kind map entry is skipped and neither the partition
// groups nor the FOR streams are walked at all.
static void test_zero_match_skips_decode() {
    int k = 5;
    uint32_t ts = table_size(k);

    // sid 1 -> short_occ1, sid 3 -> short_occ_ge2, sid 7 -> partition.
    std::vector<KpxWriter::PostingEntry> entries = {
        {1, 10},
        {3, 20}, {3, 21}, {3, 22},
    };
    for (uint32_t j = 0; j < 12; j++) entries.push_back({7, 40 + j});  // 12 > 8

    {
        KpxWriter writer(k, /*freq_threshold_part=*/8);
        for (uint32_t i = 0; i < ts; i++) {
            if (i == 31) writer.add_posting_list(i, entries);
            else         writer.add_posting_list(i, {});
        }
        CHECK(writer.write(TEST_FILE));
    }

    {
        KpxReader reader;
        CHECK(reader.open(TEST_FILE));

        const std::vector<uint32_t> kix_decoded = {1, 3, 7};
        const std::vector<uint32_t> candidates  = {2, 4, 5};

        pfd::PosDecodeScratch scratch;
        const uint64_t off = reader.pos_offset(31);
        CHECK(pfd::open_stream_kpx_for_candidates(
            reader.posting_file() + off,
            static_cast<size_t>(reader.posting_file_size() - off),
            kix_decoded.data(), kix_decoded.size(),
            candidates.data(), candidates.size(),
            scratch));

        CHECK_EQ(scratch.out_candidate_idx.size(), 0u);
        CHECK_EQ(scratch.out_positions.size(), 0u);
        CHECK_EQ(scratch.skipped_kind_entries, kix_decoded.size());
        // Nothing past the kind map was read, so nothing was walked past.
        CHECK_EQ(scratch.skipped_partition_groups, 0u);
        CHECK_EQ(scratch.skipped_for_blocks, 0u);
        reader.close();
    }
    std::remove(TEST_FILE);
}

// FOR blocks holding no selected element are walked past through their
// header.  300 short_occ1 entries = 2 full blocks + a 44-element tail; both
// candidates sit in the first block, so the second block and the tail are
// skipped.
static void test_for_block_skip() {
    int k = 5;
    uint32_t ts = table_size(k);

    const uint32_t n = 300;
    std::vector<KpxWriter::PostingEntry> entries;
    std::vector<uint32_t> kix_decoded;
    entries.reserve(n);
    kix_decoded.reserve(n);
    for (uint32_t sid = 0; sid < n; sid++) {
        entries.push_back({sid, 1000u + sid * 3u});
        kix_decoded.push_back(sid);
    }

    {
        KpxWriter writer(k);
        for (uint32_t i = 0; i < ts; i++) {
            if (i == 13) writer.add_posting_list(i, entries);
            else         writer.add_posting_list(i, {});
        }
        CHECK(writer.write(TEST_FILE));
    }

    {
        KpxReader reader;
        CHECK(reader.open(TEST_FILE));

        const std::vector<uint32_t> candidates = {0, 5};  // both in block 0

        pfd::PosDecodeScratch scratch;
        const uint64_t off = reader.pos_offset(13);
        CHECK(pfd::open_stream_kpx_for_candidates(
            reader.posting_file() + off,
            static_cast<size_t>(reader.posting_file_size() - off),
            kix_decoded.data(), kix_decoded.size(),
            candidates.data(), candidates.size(),
            scratch));

        const std::vector<uint32_t> want_idx = {0, 1};
        const std::vector<uint32_t> want_off = {0, 1, 2};
        CHECK(scratch.out_candidate_idx == want_idx);
        CHECK(scratch.out_offsets == want_off);
        CHECK_EQ(scratch.out_positions.size(), 2u);
        CHECK_EQ(scratch.out_positions[0], 1000u);
        CHECK_EQ(scratch.out_positions[1], 1015u);
        // Block 1 and the 44-element tail carry no selected element.
        CHECK_EQ(scratch.skipped_for_blocks, 2u);
        CHECK_EQ(scratch.skipped_partition_groups, 0u);
        reader.close();
    }
    std::remove(TEST_FILE);
}

// Partition groups no candidate asked for have their position stream walked
// past; only the group header is read to find the next one.
static void test_partition_group_skip() {
    int k = 5;
    uint32_t ts = table_size(k);

    const std::vector<uint32_t> sids = {2, 4, 6, 8};
    std::vector<KpxWriter::PostingEntry> entries;
    for (uint32_t sid : sids) {
        for (uint32_t j = 0; j < 12; j++) {   // 12 > 8 -> partition
            entries.push_back({sid, sid * 1000u + j});
        }
    }

    {
        KpxWriter writer(k, /*freq_threshold_part=*/8);
        for (uint32_t i = 0; i < ts; i++) {
            if (i == 27) writer.add_posting_list(i, entries);
            else         writer.add_posting_list(i, {});
        }
        CHECK(writer.write(TEST_FILE));
    }

    {
        KpxReader reader;
        CHECK(reader.open(TEST_FILE));

        const std::vector<uint32_t> candidates = {4};  // partition rank 1 only

        pfd::PosDecodeScratch scratch;
        const uint64_t off = reader.pos_offset(27);
        CHECK(pfd::open_stream_kpx_for_candidates(
            reader.posting_file() + off,
            static_cast<size_t>(reader.posting_file_size() - off),
            sids.data(), sids.size(),
            candidates.data(), candidates.size(),
            scratch));

        const std::vector<uint32_t> want_idx = {0};
        const std::vector<uint32_t> want_off = {0, 12};
        CHECK(scratch.out_candidate_idx == want_idx);
        CHECK(scratch.out_offsets == want_off);
        for (uint32_t j = 0; j < 12; j++) {
            CHECK_EQ(scratch.out_positions[j], 4000u + j);
        }
        CHECK_EQ(scratch.skipped_partition_groups, 3u);
        // Each skipped group is a 12-element stream: no full block, one tail.
        CHECK_EQ(scratch.skipped_for_blocks, 3u);
        reader.close();
    }
    std::remove(TEST_FILE);
}

// A sparse candidate set leaves most of the kind map unvisited: only the
// entries the two arrays actually meet at are resolved one at a time, the
// rest are jumped over and their ranks recovered in bulk.
static void test_kind_entry_skip() {
    int k = 5;
    uint32_t ts = table_size(k);

    const uint32_t n = 200;
    std::vector<KpxWriter::PostingEntry> entries;
    std::vector<uint32_t> kix_decoded;
    entries.reserve(n);
    kix_decoded.reserve(n);
    for (uint32_t sid = 0; sid < n; sid++) {
        entries.push_back({sid, 7000u + sid * 5u});
        kix_decoded.push_back(sid);
    }

    {
        KpxWriter writer(k);
        for (uint32_t i = 0; i < ts; i++) {
            if (i == 41) writer.add_posting_list(i, entries);
            else         writer.add_posting_list(i, {});
        }
        CHECK(writer.write(TEST_FILE));
    }

    {
        KpxReader reader;
        CHECK(reader.open(TEST_FILE));

        const std::vector<uint32_t> candidates = {0, 100, 199};

        pfd::PosDecodeScratch scratch;
        const uint64_t off = reader.pos_offset(41);
        CHECK(pfd::open_stream_kpx_for_candidates(
            reader.posting_file() + off,
            static_cast<size_t>(reader.posting_file_size() - off),
            kix_decoded.data(), kix_decoded.size(),
            candidates.data(), candidates.size(),
            scratch));

        const std::vector<uint32_t> want_idx = {0, 1, 2};
        const std::vector<uint32_t> want_off = {0, 1, 2, 3};
        CHECK(scratch.out_candidate_idx == want_idx);
        CHECK(scratch.out_offsets == want_off);
        CHECK_EQ(scratch.out_positions[0], 7000u);
        CHECK_EQ(scratch.out_positions[1], 7000u + 100u * 5u);
        CHECK_EQ(scratch.out_positions[2], 7000u + 199u * 5u);
        CHECK_EQ(scratch.skipped_kind_entries, n - candidates.size());
        reader.close();
    }
    std::remove(TEST_FILE);
}

static void test_multiple_kmers() {
    int k = 5;
    uint32_t ts = table_size(k);

    std::vector<KpxWriter::PostingEntry> entries_a = {
        {0, 100}, {0, 200}, {1, 50},
    };
    std::vector<uint32_t> cand_a = {0, 1};

    std::vector<KpxWriter::PostingEntry> entries_b = {
        {2, 0}, {2, 1000},
    };
    std::vector<uint32_t> cand_b = {2};

    {
        KpxWriter writer(k);
        for (uint32_t i = 0; i < ts; i++) {
            if (i == 0) {
                writer.add_posting_list(i, entries_a);
            } else if (i == 10) {
                writer.add_posting_list(i, entries_b);
            } else {
                writer.add_posting_list(i, {});
            }
        }
        CHECK(writer.write(TEST_FILE));
    }

    {
        KpxReader reader;
        CHECK(reader.open(TEST_FILE));
        CHECK_EQ(reader.total_position_count(), 5u);

        auto pos_a = decode_pos_postings(
            reader.posting_file() + reader.pos_offset(0),
            reader.posting_file() + reader.posting_file_size(),
            cand_a);
        CHECK_EQ(pos_a.size(), 3u);
        CHECK_EQ(pos_a[0], 100u);
        CHECK_EQ(pos_a[1], 200u);
        CHECK_EQ(pos_a[2], 50u);

        auto pos_b = decode_pos_postings(
            reader.posting_file() + reader.pos_offset(10),
            reader.posting_file() + reader.posting_file_size(),
            cand_b);
        CHECK_EQ(pos_b.size(), 2u);
        CHECK_EQ(pos_b[0], 0u);
        CHECK_EQ(pos_b[1], 1000u);

        reader.close();
    }

    std::remove(TEST_FILE);
}

// Partition / short bucket boundary across thresholds.  A 20-occurrence
// run on seq_id=10 plus 1 occurrence each on {11, 12, 13} (occ=1 →
// short_occ1 sub-bucket).
//   threshold = 8   -> seq_id=10 partitioned, others short_occ1
//   threshold = 16  -> seq_id=10 partitioned (20 > 16), others short_occ1
//   threshold = 100 -> all clusters short_occ1, no partition groups
static void test_partition_group_threshold() {
    int k = 5;
    uint32_t ts = table_size(k);

    std::vector<KpxWriter::PostingEntry> entries;
    for (uint32_t pos = 0; pos < 20; pos++) {
        entries.push_back({10, 1000 + pos});
    }
    for (uint32_t sid : {11u, 12u, 13u}) {
        entries.push_back({sid, 7});
    }
    const uint32_t expected_count = static_cast<uint32_t>(entries.size());
    std::vector<uint32_t> candidates = {10, 11, 12, 13};

    auto run_threshold = [&](uint32_t threshold) {
        KpxWriter writer(k, threshold);
        for (uint32_t i = 0; i < ts; i++) {
            if (i == 33) writer.add_posting_list(i, entries);
            else        writer.add_posting_list(i, {});
        }
        CHECK(writer.write(TEST_FILE));

        KpxReader reader;
        CHECK(reader.open(TEST_FILE));
        CHECK_EQ(reader.total_position_count(), expected_count);

        auto positions = decode_pos_postings(
            reader.posting_file() + reader.pos_offset(33),
            reader.posting_file() + reader.posting_file_size(),
            candidates);

        CHECK_EQ(positions.size(), static_cast<size_t>(expected_count));
        for (size_t j = 0; j < entries.size(); j++) {
            CHECK_EQ(positions[j], entries[j].pos);
        }

        reader.close();
        std::remove(TEST_FILE);
    };

    run_threshold(8);    // seq_id=10 cluster (20 > 8) → partition
    run_threshold(16);   // seq_id=10 cluster (20 > 16) → partition
    run_threshold(100);  // nothing exceeds 100 → all short
}

int main() {
    test_single_seq();
    test_all_occ1();
    test_mixed_short_buckets();
    test_partition_only();
    test_partition_tail_sizes();
    test_partition_above_255();
    test_short1_full_blocks();
    test_empty_posting();
    test_all_candidates_miss();
    test_candidate_subset_csr();
    test_zero_match_skips_decode();
    test_for_block_skip();
    test_partition_group_skip();
    test_kind_entry_skip();
    test_multiple_kmers();
    test_partition_group_threshold();
    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
