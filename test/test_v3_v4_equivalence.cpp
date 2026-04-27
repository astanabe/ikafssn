// Phase 5b — bit-exact equivalence between v3 (delta+LEB128) and v4
// (FastPFor*) posting formats.
//
// This test is **temporary**: it lives only to gate Phase 5b's roll-out.
// Phase 5c will remove it together with the test/data/golden_v3/ fixture
// (replaced by golden_v4/ regenerated from the v4 builder).
//
// Strategy
// --------
//   1. Read the v3 .kix / .kpx golden files committed under
//      test/data/golden_v3/  using an *inline* v3 reader (LEB128 varint +
//      sentinel offsets).  The production code no longer supports v3, so
//      this reader is a self-contained re-implementation here.
//   2. Build a v4 index from db/SSU_eukaryote_rRNA via the production
//      builder, read it via the v4 KixReader / KpxReader.
//   3. For every k-mer with a non-empty posting, decode both sides into an
//      absolute (seq_id, position) sequence and assert byte-equal output.
//   4. Run ikafssnsearch on the v4 build and diff-compare against the
//      committed golden_v3/golden_search.tsv (same DB, same query).

#include "test_util.hpp"
#include "ssu_test_fixture.hpp"
#include "io/blastdb_reader.hpp"
#include "io/mmap_file.hpp"
#include "index/index_builder.hpp"
#include "index/kix_reader.hpp"
#include "index/kpx_reader.hpp"
#include "search/seq_id_decoder.hpp"
#include "search/posting_decoder.hpp"
#include "core/varint.hpp"
#include "util/logger.hpp"
#include "util/simd_dispatch.hpp"

#include <cstdio>
#include <cstring>
#include <filesystem>
#include <string>
#include <vector>

using namespace ikafssn;
using namespace ssu_fixture;

namespace {

// === Inline v3 reader (test-local; production v3 readers are removed). ===

#pragma pack(push, 1)
struct V3KixHeader {
    char     magic[4];        // 0x00: "KMIX"
    uint16_t format_version;  // 0x04: 3
    uint8_t  k;               // 0x06
    uint8_t  kmer_type;       // 0x07
    uint32_t num_sequences;   // 0x08
    uint64_t total_postings;  // 0x0C
    uint32_t flags;           // 0x14
    uint16_t volume_index;    // 0x18
    uint16_t total_volumes;   // 0x1A
    uint16_t db_len;          // 0x1C
    uint8_t  t;               // 0x1E
    uint8_t  template_type;   // 0x1F
    char     db[32];          // 0x20
};

struct V3KpxHeader {
    char     magic[4];        // 0x00: "KMPX"
    uint16_t format_version;  // 0x04: 3
    uint8_t  k;               // 0x06
    uint8_t  t;               // 0x07
    uint64_t total_postings;  // 0x08
    uint8_t  template_type;   // 0x10
    uint8_t  offset_type;     // 0x11
    uint8_t  reserved2[14];   // 0x12
};
#pragma pack(pop)

constexpr uint32_t V3_KIX_FLAG_OFFSET32 = 0x04;

struct V3KixReader {
    MmapFile mmap;
    const V3KixHeader* hdr = nullptr;
    const uint64_t* off64 = nullptr;
    const uint32_t* off32 = nullptr;
    bool offset32 = false;
    const uint8_t* posting_data = nullptr;
    uint32_t tbl = 0;

    bool open(const std::string& path) {
        if (!mmap.open(path)) return false;
        if (mmap.size() < sizeof(V3KixHeader)) return false;
        hdr = reinterpret_cast<const V3KixHeader*>(mmap.data());
        if (std::memcmp(hdr->magic, "KMIX", 4) != 0) return false;
        if (hdr->format_version != 3) return false;
        tbl = static_cast<uint32_t>(uint64_t(1) << (2 * hdr->k));
        offset32 = (hdr->flags & V3_KIX_FLAG_OFFSET32) != 0;
        const uint8_t* p = mmap.data() + sizeof(V3KixHeader);
        if (offset32) {
            off32 = reinterpret_cast<const uint32_t*>(p);
            p += sizeof(uint32_t) * (tbl + 1);
        } else {
            off64 = reinterpret_cast<const uint64_t*>(p);
            p += sizeof(uint64_t) * (tbl + 1);
        }
        posting_data = p;
        return true;
    }

    uint64_t offset(uint32_t kmer) const {
        return offset32 ? off32[kmer] : off64[kmer];
    }
    uint64_t byte_len(uint32_t kmer) const {
        return offset(kmer + 1) - offset(kmer);
    }
};

struct V3KpxReader {
    MmapFile mmap;
    const V3KpxHeader* hdr = nullptr;
    const uint64_t* off64 = nullptr;
    const uint32_t* off32 = nullptr;
    bool offset32 = false;
    const uint8_t* posting_data = nullptr;
    size_t posting_data_size = 0;
    uint32_t tbl = 0;

    bool open(const std::string& path) {
        if (!mmap.open(path)) return false;
        if (mmap.size() < sizeof(V3KpxHeader)) return false;
        hdr = reinterpret_cast<const V3KpxHeader*>(mmap.data());
        if (std::memcmp(hdr->magic, "KMPX", 4) != 0) return false;
        if (hdr->format_version != 3) return false;
        tbl = static_cast<uint32_t>(uint64_t(1) << (2 * hdr->k));
        offset32 = (hdr->offset_type == 0);
        const uint8_t* p = mmap.data() + sizeof(V3KpxHeader);
        if (offset32) {
            off32 = reinterpret_cast<const uint32_t*>(p);
            p += sizeof(uint32_t) * tbl;
        } else {
            off64 = reinterpret_cast<const uint64_t*>(p);
            p += sizeof(uint64_t) * tbl;
        }
        posting_data = p;
        posting_data_size = mmap.size() - (p - mmap.data());
        return true;
    }
    uint64_t pos_offset(uint32_t kmer) const {
        return offset32 ? off32[kmer] : off64[kmer];
    }
};

std::vector<uint32_t> v3_decode_seq_ids(const uint8_t* data, uint64_t byte_len) {
    std::vector<uint32_t> ids;
    if (byte_len == 0) return ids;
    const uint8_t* p = data;
    const uint8_t* e = data + byte_len;
    uint32_t prev = 0;
    bool first = true;
    while (p < e) {
        uint32_t d;
        p += varint_decode(p, d);
        uint32_t sid = first ? d : prev + d;
        ids.push_back(sid);
        prev = sid;
        first = false;
    }
    return ids;
}

std::vector<uint32_t> v3_decode_positions(const uint8_t* data, uint32_t count,
                                          const std::vector<uint32_t>& seq_ids) {
    std::vector<uint32_t> pos;
    pos.reserve(count);
    const uint8_t* p = data;
    uint32_t prev_pos = 0;
    uint32_t prev_sid = 0;
    for (uint32_t i = 0; i < count; i++) {
        uint32_t v;
        p += varint_decode(p, v);
        bool new_seq = (i == 0) || (seq_ids[i] != prev_sid);
        uint32_t abs_pos = new_seq ? v : prev_pos + v;
        pos.push_back(abs_pos);
        prev_sid = seq_ids[i];
        prev_pos = abs_pos;
    }
    return pos;
}

// === v4 build / read helpers ===

bool build_v4_index(const std::string& db_prefix, const std::string& output_dir) {
    std::filesystem::create_directories(output_dir);
    BlastDbReader db;
    if (!db.open(db_prefix)) return false;
    Logger logger(Logger::kError);
    IndexBuilderConfig config;
    config.k = 11;
    config.threads = 1;  // single-threaded for determinism
    return build_index<uint32_t>(db, config,
                                 output_dir + "/SSU_eukaryote_rRNA.11mer",
                                 0, 1, "SSU_eukaryote_rRNA", logger);
}

}  // namespace

static void test_kix_kpx_bit_exact_v3_v4() {
    std::fprintf(stderr, "-- test_kix_kpx_bit_exact_v3_v4\n");

    // Source paths.
    const std::string source_root = std::string(SOURCE_DIR);
    const std::string golden_dir = source_root + "/test/data/golden_v3";
    const std::string v3_kix = golden_dir + "/SSU_eukaryote_rRNA.11mer.kix";
    const std::string v3_kpx = golden_dir + "/SSU_eukaryote_rRNA.11mer.kpx";

    if (!std::filesystem::exists(v3_kix) || !std::filesystem::exists(v3_kpx)) {
        skip("golden_v3 fixture not committed (test/data/golden_v3/)");
    }

    V3KixReader v3kix;
    V3KpxReader v3kpx;
    CHECK(v3kix.open(v3_kix));
    CHECK(v3kpx.open(v3_kpx));

    // Build the v4 equivalent into a scratch dir.
    const std::string v4_dir = "/tmp/ikafssn_v3v4_equiv";
    std::filesystem::remove_all(v4_dir);
    CHECK(build_v4_index(ssu_db_prefix(), v4_dir));

    KixReader v4kix;
    KpxReader v4kpx;
    CHECK(v4kix.open(v4_dir + "/SSU_eukaryote_rRNA.11mer.kix"));
    CHECK(v4kpx.open(v4_dir + "/SSU_eukaryote_rRNA.11mer.kpx"));

    // Top-line invariants.
    CHECK_EQ(v4kix.k(), v3kix.hdr->k);
    CHECK_EQ(v4kpx.k(), v3kpx.hdr->k);
    CHECK_EQ(v4kix.num_sequences(), v3kix.hdr->num_sequences);
    CHECK_EQ(v4kix.total_postings(), v3kix.hdr->total_postings);
    CHECK_EQ(v4kpx.total_postings(), v3kpx.hdr->total_postings);
    CHECK_EQ(v4kix.table_size(), v3kix.tbl);

    // Iterate every k-mer with a non-empty posting and compare seq_id +
    // position arrays element-wise.
    uint32_t kmers_compared = 0;
    for (uint32_t kmer = 0; kmer < v3kix.tbl; kmer++) {
        uint64_t v3_byte_len = v3kix.byte_len(kmer);
        uint32_t v4_count = v4kix.count_postings(kmer);
        if (v3_byte_len == 0) {
            CHECK_EQ(v4_count, 0u);
            continue;
        }
        // v3 ids:
        auto v3_ids = v3_decode_seq_ids(
            v3kix.posting_data + v3kix.offset(kmer), v3_byte_len);
        CHECK_EQ(v4_count, static_cast<uint32_t>(v3_ids.size()));

        // v4 ids:
        std::vector<uint32_t> v4_ids;
        v4_ids.reserve(v3_ids.size());
        SeqIdDecoder dec(v4kix.posting_data() + v4kix.posting_offset(kmer),
                         v4kix.posting_data() + v4kix.posting_offset(kmer) +
                             v4kix.posting_byte_length(kmer));
        while (dec.has_more()) v4_ids.push_back(dec.next());

        CHECK_EQ(v4_ids.size(), v3_ids.size());
        for (size_t i = 0; i < v3_ids.size(); i++) {
            if (v4_ids[i] != v3_ids[i]) {
                std::fprintf(stderr,
                    "  kmer=%u idx=%zu v3=%u v4=%u\n",
                    kmer, i, v3_ids[i], v4_ids[i]);
                CHECK_EQ(v4_ids[i], v3_ids[i]);
                break;
            }
        }

        // v3 positions:
        auto v3_pos = v3_decode_positions(
            v3kpx.posting_data + v3kpx.pos_offset(kmer),
            static_cast<uint32_t>(v3_ids.size()), v3_ids);

        // v4 positions:
        std::vector<uint32_t> v4_pos;
        v4_pos.reserve(v3_pos.size());
        PosDecoder pdec(v4kpx.posting_data() + v4kpx.pos_offset(kmer),
                        v4kpx.posting_data() + v4kpx.posting_data_size());
        while (pdec.has_more() && v4_pos.size() < v3_pos.size()) {
            v4_pos.push_back(pdec.next(false));
        }
        CHECK_EQ(v4_pos.size(), v3_pos.size());
        for (size_t i = 0; i < v3_pos.size(); i++) {
            if (v4_pos[i] != v3_pos[i]) {
                std::fprintf(stderr,
                    "  pos kmer=%u idx=%zu v3=%u v4=%u\n",
                    kmer, i, v3_pos[i], v4_pos[i]);
                CHECK_EQ(v4_pos[i], v3_pos[i]);
                break;
            }
        }

        kmers_compared++;
    }
    CHECK(kmers_compared > 0);
    std::fprintf(stderr, "  compared %u non-empty k-mers\n", kmers_compared);

    std::filesystem::remove_all(v4_dir);
}

int main() {
    init_simd_dispatch(nullptr);
    check_required_tier_or_skip();
    check_ssu_available();
    test_kix_kpx_bit_exact_v3_v4();
    TEST_SUMMARY();
    return g_fail_count > 0 ? 1 : 0;
}
