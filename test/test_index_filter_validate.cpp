#include "test_util.hpp"
#include "kpx_writer.hpp"
#include "index/kpx_reader.hpp"
#include "index/kpx_format.hpp"
#include "kix_writer.hpp"
#include "index/index_filter.hpp"
#include "core/config.hpp"
#include "util/logger.hpp"

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

using namespace ikafssn;

namespace {

// Build a tiny synthetic (.kix, .kpx) volume pair at `prefix`.{kix,kpx}.
// The k=5 dictionary has two non-empty k-mers, one large enough to
// exercise the partition path and one short_occ_ge2 cluster.  Returns
// the table size.
uint32_t build_synth_volume(const std::string& prefix) {
    const int k = 5;
    const uint32_t ts = table_size(k);

    // .kix: two k-mers with some seq_ids.
    KixWriter kix(k, /*kmer_type=*/0);
    kix.set_num_sequences(64);
    for (uint32_t i = 0; i < ts; i++) {
        if (i == 7) {
            kix.add_posting_list(i, {1, 2, 3, 5, 7, 11, 13, 17, 19, 23});
        } else if (i == 42) {
            kix.add_posting_list(i, {3, 4});
        } else {
            kix.add_posting_list(i, {});
        }
    }
    CHECK(kix.write(prefix + ".kix"));

    // .kpx: matching distinct sids.  freq_threshold_part = 8 default;
    // the kmer=7 cluster (sid=1, 100 occurrences) lands as a partition
    // group; sid=2 lands as short_occ1 (single position); the rest go
    // into short_occ_ge2.
    KpxWriter kpx(k);
    for (uint32_t i = 0; i < ts; i++) {
        if (i == 7) {
            std::vector<KpxWriter::PostingEntry> entries;
            // sid=1: partition group (>8 occurrences)
            for (uint32_t pos = 0; pos < 50; pos++) {
                entries.push_back({1, 1000u + pos});
            }
            // sid=2: short_occ1
            entries.push_back({2, 555});
            // sid=3..23: short_occ_ge2 with 2-3 positions each
            for (uint32_t sid : {3u, 5u, 7u, 11u, 13u, 17u, 19u, 23u}) {
                entries.push_back({sid, 100u + sid});
                entries.push_back({sid, 200u + sid});
            }
            kpx.add_posting_list(i, entries);
        } else if (i == 42) {
            std::vector<KpxWriter::PostingEntry> entries = {
                {3, 17}, {4, 19}
            };
            kpx.add_posting_list(i, entries);
        } else {
            kpx.add_posting_list(i, {});
        }
    }
    CHECK(kpx.write(prefix + ".kpx"));

    return ts;
}

void test_validate_clean_volume() {
    std::string prefix = test_tmpdir("/tmp/test_validate_clean") + "/vol";
    std::filesystem::remove_all(test_tmpdir("/tmp/test_validate_clean"));
    std::filesystem::create_directories(test_tmpdir("/tmp/test_validate_clean"));
    build_synth_volume(prefix);

    Logger logger(Logger::kError);
    uint64_t total_pos = 0;
    bool ok = validate_volume(prefix + ".kix", prefix + ".kpx", &total_pos, logger);
    CHECK(ok);
    // sid=1 partition (50) + sid=2 short_occ1 (1) + 8 * 2 short_occ_ge2 (16)
    //   for kmer=7 -> 67 positions
    // kmer=42: 2 positions
    // total = 69
    CHECK_EQ(total_pos, 69u);

    std::filesystem::remove_all(test_tmpdir("/tmp/test_validate_clean"));
}

// Inject a single-bit flip into one of the kind-map bytes of a known
// k-mer's .kpx posting list and confirm the validator catches it.
void test_validate_detects_kind_map_corruption() {
    std::string prefix = test_tmpdir("/tmp/test_validate_corrupt") + "/vol";
    std::filesystem::remove_all(test_tmpdir("/tmp/test_validate_corrupt"));
    std::filesystem::create_directories(test_tmpdir("/tmp/test_validate_corrupt"));
    build_synth_volume(prefix);

    // Sanity: the clean file passes.
    {
        Logger logger(Logger::kError);
        CHECK(validate_volume(prefix + ".kix", prefix + ".kpx", nullptr, logger));
    }

    // Inject a corruption: flip the high-order bit of the byte at the
    // start of kmer=7's .kpx posting list.  The validator opens the
    // file via mmap so we have to perform the flip on disk and reopen.
    {
        std::size_t abs_off = 0;
        {
            KpxReader reader;
            CHECK(reader.open(prefix + ".kpx"));
            const uint64_t kmer7_off = reader.pos_offset(7);
            const std::ptrdiff_t posting_start_in_file =
                reader.posting_file() - reinterpret_cast<const uint8_t*>(&reader.header());
            abs_off = static_cast<std::size_t>(posting_start_in_file) + kmer7_off;
            reader.close();
        }

        std::fstream fs(prefix + ".kpx",
                        std::ios::in | std::ios::out | std::ios::binary);
        CHECK(fs.is_open());
        fs.seekg(static_cast<std::streamoff>(abs_off));
        char b;
        fs.read(&b, 1);
        b ^= static_cast<char>(0x80);
        fs.seekp(static_cast<std::streamoff>(abs_off));
        fs.write(&b, 1);
        fs.close();
    }

    // The validator should now reject the volume.  The expected error
    // message is intentionally allowed through to stderr.
    {
        Logger logger(Logger::kError);
        bool ok = validate_volume(prefix + ".kix", prefix + ".kpx", nullptr, logger);
        CHECK(!ok);
    }

    std::filesystem::remove_all(test_tmpdir("/tmp/test_validate_corrupt"));
}

// Validate a mode-1 build (.kix only, no .kpx) by passing an empty
// kpx_path.  Should pass the structural check unconditionally because
// the .kix-only walk has no posting list bodies to verify.
void test_validate_mode1_kix_only() {
    std::string prefix = test_tmpdir("/tmp/test_validate_mode1") + "/vol";
    std::filesystem::remove_all(test_tmpdir("/tmp/test_validate_mode1"));
    std::filesystem::create_directories(test_tmpdir("/tmp/test_validate_mode1"));
    build_synth_volume(prefix);

    Logger logger(Logger::kError);
    CHECK(validate_volume(prefix + ".kix", "", nullptr, logger));

    std::filesystem::remove_all(test_tmpdir("/tmp/test_validate_mode1"));
}

} // anonymous namespace

int main() {
    test_validate_clean_volume();
    test_validate_detects_kind_map_corruption();
    test_validate_mode1_kix_only();
    TEST_SUMMARY();
    return g_fail_count == 0 ? 0 : 1;
}
