// Compression analyzer for .kix and .kpx posting lists.
//
// Walks the binary directly (no decode into vectors) and reports:
//
//   - Bytes attributed to: top headers, kind map bytes, partition group
//     headers, short_occ1 / short_occ_ge2 sub-bucket bytes, FOR block
//     headers, bitpacked block bodies, packed-bit tail bodies, .kix
//     posting list body bytes.
//   - Per-block bit-width 'b' histogram, split by partition / short_occ1
//     / short_occ_ge2.
//   - Per-(partition / short_occ1 / short_occ_ge2 / .kix) bytes-per-element.
//   - Element counts.
//
// Usage:
//   index_compress_stats --kix <path>            # .kix only
//   index_compress_stats --kpx <path>            # .kpx only
//   index_compress_stats --kix <p> --kpx <p>     # both
//
// The tool does not parse the body words of `.kix` posting lists; it only
// sums them so we can measure bytes-per-distinct-seq_id.  Intra-block 'b'
// for .kpx is read directly from the on-disk byte stream.

#include "index/kix_reader.hpp"
#include "index/kpx_reader.hpp"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <array>

using namespace ikafssn;

namespace {

constexpr int kBlockSize = 128;

struct ForStreamStats {
    uint64_t blocks = 0;
    uint64_t block_header_bytes = 0;   // 8 B per full block
    uint64_t block_body_bytes = 0;     // 16*b per full block
    uint64_t tail_count_bytes = 0;     // 1 B per stream
    uint64_t tail_min_bytes = 0;       // 4 B when tail_count > 0
    uint64_t tail_b_bytes = 0;         // 1 B when tail_count > 0
    uint64_t tail_body_bytes = 0;      // ceil(tail_count*tail_b/8)
    uint64_t tail_elems = 0;
    std::array<uint64_t, 33> b_hist{};
    std::array<uint64_t, 33> tail_b_hist{};
};

inline uint64_t for_stream_total(const ForStreamStats& s) {
    return s.block_header_bytes + s.block_body_bytes +
           s.tail_count_bytes + s.tail_min_bytes +
           s.tail_b_bytes + s.tail_body_bytes;
}

struct KpxStats {
    uint64_t n_kmers_nonempty = 0;
    uint64_t total_position_count = 0;
    uint64_t total_distinct_count = 0;

    // Reserved field, always 0 (no top-level posting list header).
    uint64_t bytes_top_header = 0;
    // 2-bit kind map bytes per posting list.
    uint64_t bytes_kind_map = 0;

    // Partition groups.
    uint64_t partition_groups = 0;
    uint64_t partition_positions = 0;
    uint64_t partition_header_bytes = 0;       // 4 B per group (occ_count)
    ForStreamStats partition_stream;
    uint64_t partition_occ_min = UINT64_MAX, partition_occ_max = 0;
    uint64_t partition_occ_sum = 0;

    // Short occ=1 sub-bucket.
    uint64_t short1_seqs = 0;
    uint64_t short1_positions = 0;
    ForStreamStats short1_stream;

    // Short occ>=2 sub-bucket.
    uint64_t short2_seqs = 0;
    uint64_t short2_positions = 0;
    uint64_t short2_occ_bytes = 0;             // u8 per seq
    ForStreamStats short2_stream;
};

struct KixStats {
    uint64_t n_kmers_nonempty = 0;
    uint64_t total_distinct_count = 0;
    // Per posting list: 8 B header (count + body_words) + body_words*4.
    uint64_t bytes_posting_list_header = 0;
    uint64_t bytes_posting_list_body = 0;
};

// Walk a FOR block stream of `count` elements, accumulating stats.
// Returns advanced bytes; -1 on corruption.
ssize_t walk_for_stream(const uint8_t* p, const uint8_t* end, uint32_t count,
                        ForStreamStats& s) {
    const uint8_t* p0 = p;
    const uint32_t num_blocks = count / kBlockSize;
    const uint32_t tail_count = count % kBlockSize;
    for (uint32_t bi = 0; bi < num_blocks; bi++) {
        if (size_t(end - p) < 8) return -1;
        uint8_t b = p[4];
        if (b > 32) return -1;
        s.b_hist[b]++;
        const std::size_t body_bytes = std::size_t(16) * b;
        if (size_t(end - p) < 8 + body_bytes) return -1;
        s.block_header_bytes += 8;
        s.block_body_bytes += body_bytes;
        p += 8 + body_bytes;
        s.blocks++;
    }
    if (p >= end) return -1;
    uint8_t got_tail = *p++;
    s.tail_count_bytes += 1;
    if (got_tail != tail_count) return -1;
    if (tail_count > 0) {
        if (size_t(end - p) < 5) return -1;
        p += 4;  // tail_min
        s.tail_min_bytes += 4;
        uint8_t tail_b = *p++;
        if (tail_b > 32) return -1;
        s.tail_b_bytes += 1;
        s.tail_b_hist[tail_b]++;
        const std::size_t body_bits  = std::size_t(tail_count) * tail_b;
        const std::size_t body_bytes = (body_bits + 7) / 8;
        if (size_t(end - p) < body_bytes) return -1;
        p += body_bytes;
        s.tail_body_bytes += body_bytes;
        s.tail_elems += tail_count;
    }
    return p - p0;
}

// distinct_count comes from .kix; the body starts directly at the
// 2-bit kind map.
bool walk_kpx_posting(const uint8_t* p, const uint8_t* end, uint32_t kix_count,
                      KpxStats& st) {
    if (kix_count == 0) {
        // Empty .kix posting list aliases to the next non-empty
        // k-mer's slice; the .kpx writer emits zero bytes for it.
        return true;
    }

    st.n_kmers_nonempty++;
    st.total_distinct_count += kix_count;

    // 2-bit kind map.
    const std::size_t kind_map_bytes = (std::size_t(kix_count) * 2 + 7) / 8;
    if (size_t(end - p) < kind_map_bytes) return false;
    const uint8_t* kind_map = p;
    p += kind_map_bytes;
    st.bytes_kind_map += kind_map_bytes;

    // Per-kind counts via kind-map popcount.
    uint32_t partition_count = 0, short1_count = 0, short2_count = 0;
    for (uint32_t k = 0; k < kix_count; k++) {
        uint8_t kind = static_cast<uint8_t>(
            (kind_map[k >> 2] >> ((k & 3) * 2)) & 0x03);
        switch (kind) {
            case 0: short1_count++;    break;
            case 1: short2_count++;    break;
            case 2: partition_count++; break;
            default: return false;  // 11 reserved
        }
    }

    // Partition groups: each [u32 occ_count][FOR stream].
    for (uint32_t g = 0; g < partition_count; g++) {
        if (size_t(end - p) < 4) return false;
        uint32_t gcnt;
        std::memcpy(&gcnt, p, 4); p += 4;
        st.partition_header_bytes += 4;
        st.partition_groups++;
        st.partition_positions += gcnt;
        if (gcnt < st.partition_occ_min) st.partition_occ_min = gcnt;
        if (gcnt > st.partition_occ_max) st.partition_occ_max = gcnt;
        st.partition_occ_sum += gcnt;
        ssize_t adv = walk_for_stream(p, end, gcnt, st.partition_stream);
        if (adv < 0) return false;
        p += adv;
    }

    // short_occ1 sub-bucket.
    if (short1_count > 0) {
        st.short1_seqs += short1_count;
        st.short1_positions += short1_count;  // 1 position per cluster
        ssize_t adv = walk_for_stream(p, end, short1_count, st.short1_stream);
        if (adv < 0) return false;
        p += adv;
    }

    // short_occ_ge2 sub-bucket.
    if (short2_count > 0) {
        st.short2_seqs += short2_count;
        if (size_t(end - p) < short2_count) return false;
        uint32_t short2_position_count = 0;
        for (uint32_t i = 0; i < short2_count; i++) {
            short2_position_count += p[i];
        }
        st.short2_positions += short2_position_count;
        p += short2_count;
        st.short2_occ_bytes += short2_count;
        if (short2_position_count > 0) {
            ssize_t adv = walk_for_stream(p, end, short2_position_count, st.short2_stream);
            if (adv < 0) return false;
            p += adv;
        }
    }

    return true;
}

void print_for_stream(const char* label, const ForStreamStats& s,
                      uint64_t total_positions) {
    const uint64_t total = for_stream_total(s);
    std::printf("    [%s FOR stream] total bytes: %lu\n", label, total);
    if (total_positions > 0) {
        std::printf("      bytes/pos               : %.3f (%.2f bits/pos)\n",
                    double(total) / double(total_positions),
                    double(total) * 8.0 / double(total_positions));
    }
    std::printf("      block header bytes      : %lu (8 B/block)\n", s.block_header_bytes);
    std::printf("      block body bytes        : %lu (16*b/block)\n", s.block_body_bytes);
    std::printf("      tail count bytes        : %lu\n", s.tail_count_bytes);
    std::printf("      tail min bytes          : %lu\n", s.tail_min_bytes);
    std::printf("      tail b bytes            : %lu\n", s.tail_b_bytes);
    std::printf("      tail body bytes         : %lu\n", s.tail_body_bytes);
    std::printf("      full blocks             : %lu (covering %lu positions)\n",
                s.blocks, s.blocks * uint64_t(kBlockSize));
    std::printf("      tail positions          : %lu\n", s.tail_elems);
    if (s.blocks > 0) {
        double avg_b = 0;
        for (int b = 0; b <= 32; b++) avg_b += double(b) * double(s.b_hist[b]);
        avg_b /= double(s.blocks);
        std::printf("      avg block-b             : %.2f bits\n", avg_b);
        std::printf("      block-b histogram:\n");
        for (int b = 0; b <= 32; b++) {
            if (s.b_hist[b] > 0) {
                std::printf("        b=%2d : %lu (%.1f%%)\n", b, s.b_hist[b],
                            100.0 * double(s.b_hist[b]) / double(s.blocks));
            }
        }
    }
    if (s.tail_elems > 0) {
        std::printf("      tail-b histogram:\n");
        for (int b = 0; b <= 32; b++) {
            if (s.tail_b_hist[b] > 0) {
                std::printf("        tail_b=%2d : %lu\n", b, s.tail_b_hist[b]);
            }
        }
    }
}

void print_kpx_report(const KpxStats& st, uint64_t total_bytes) {
    std::printf("\n=== .kpx compression breakdown ===\n");
    std::printf("non-empty kmers              : %lu\n", st.n_kmers_nonempty);
    std::printf("total positions              : %lu\n", st.total_position_count);
    std::printf("total distinct seq_ids       : %lu\n", st.total_distinct_count);
    std::printf("on-disk posting file bytes        : %lu (%.3f bytes/pos, %.2f bits/pos)\n",
                total_bytes,
                st.total_position_count ? double(total_bytes) / double(st.total_position_count) : 0.0,
                st.total_position_count ? double(total_bytes) * 8.0 / double(st.total_position_count) : 0.0);

    const uint64_t sum_partition = st.partition_header_bytes + for_stream_total(st.partition_stream);
    const uint64_t sum_short1    = for_stream_total(st.short1_stream);
    const uint64_t sum_short2    = st.short2_occ_bytes + for_stream_total(st.short2_stream);
    const uint64_t top_total     = st.bytes_top_header + st.bytes_kind_map;

    std::printf("\n[top headers + kind map]      : %lu (%.2f%% of file)\n",
                top_total, 100.0 * double(top_total) / double(total_bytes));
    std::printf("  top header bytes (reserved, always 0)  : %lu\n", st.bytes_top_header);
    std::printf("  2-bit kind map bytes       : %lu\n", st.bytes_kind_map);

    std::printf("\n[partition groups] total      : %lu (%.2f%% of file)\n",
                sum_partition, 100.0 * double(sum_partition) / double(total_bytes));
    std::printf("  groups                     : %lu\n", st.partition_groups);
    std::printf("  positions                  : %lu (%.2f%% of all positions)\n",
                st.partition_positions,
                st.total_position_count ? 100.0 * double(st.partition_positions) / double(st.total_position_count) : 0.0);
    if (st.partition_groups > 0) {
        std::printf("  occ count: min=%lu max=%lu mean=%.1f\n",
                    st.partition_occ_min, st.partition_occ_max,
                    double(st.partition_occ_sum) / double(st.partition_groups));
    }
    if (st.partition_positions > 0) {
        std::printf("  bytes/pos                  : %.3f (%.2f bits/pos)\n",
                    double(sum_partition) / double(st.partition_positions),
                    double(sum_partition) * 8.0 / double(st.partition_positions));
    }
    std::printf("    group header bytes (occ) : %lu (4 B/group)\n", st.partition_header_bytes);
    print_for_stream("partition", st.partition_stream, st.partition_positions);

    std::printf("\n[short_occ1 sub-bucket] total : %lu (%.2f%% of file)\n",
                sum_short1, 100.0 * double(sum_short1) / double(total_bytes));
    std::printf("  short1 seqs                : %lu\n", st.short1_seqs);
    std::printf("  short1 positions           : %lu (%.2f%% of all positions)\n",
                st.short1_positions,
                st.total_position_count ? 100.0 * double(st.short1_positions) / double(st.total_position_count) : 0.0);
    if (st.short1_positions > 0) {
        std::printf("  bytes/pos                  : %.3f (%.2f bits/pos)\n",
                    double(sum_short1) / double(st.short1_positions),
                    double(sum_short1) * 8.0 / double(st.short1_positions));
    }
    print_for_stream("short_occ1", st.short1_stream, st.short1_positions);

    std::printf("\n[short_occ_ge2 sub-bucket]    : %lu (%.2f%% of file)\n",
                sum_short2, 100.0 * double(sum_short2) / double(total_bytes));
    std::printf("  short2 seqs                : %lu\n", st.short2_seqs);
    std::printf("  short2 positions           : %lu (%.2f%% of all positions)\n",
                st.short2_positions,
                st.total_position_count ? 100.0 * double(st.short2_positions) / double(st.total_position_count) : 0.0);
    if (st.short2_seqs > 0) {
        std::printf("  positions/short2_seq       : %.3f\n",
                    double(st.short2_positions) / double(st.short2_seqs));
    }
    if (st.short2_positions > 0) {
        std::printf("  bytes/pos                  : %.3f (%.2f bits/pos)\n",
                    double(sum_short2) / double(st.short2_positions),
                    double(sum_short2) * 8.0 / double(st.short2_positions));
    }
    std::printf("    occ_count bytes (u8)     : %lu\n", st.short2_occ_bytes);
    print_for_stream("short_occ_ge2", st.short2_stream, st.short2_positions);

    // Sanity: top + partition + short1 + short2 should equal total file.
    const uint64_t reconstructed = top_total + sum_partition + sum_short1 + sum_short2;
    if (reconstructed != total_bytes) {
        std::printf("\nNOTE: walked sum %lu vs file size %lu (delta %ld)\n",
                    reconstructed, total_bytes, long(total_bytes) - long(reconstructed));
    }
}

void print_kix_report(const KixStats& st, uint64_t total_bytes) {
    std::printf("\n=== .kix compression breakdown ===\n");
    std::printf("non-empty kmers              : %lu\n", st.n_kmers_nonempty);
    std::printf("total distinct seq_ids       : %lu\n", st.total_distinct_count);
    std::printf("on-disk posting list bytes   : %lu\n", total_bytes);
    if (st.total_distinct_count > 0) {
        std::printf("  bytes/distinct_seq_id      : %.3f (%.2f bits/sid)\n",
                    double(total_bytes) / double(st.total_distinct_count),
                    double(total_bytes) * 8.0 / double(st.total_distinct_count));
    }
    std::printf("  posting list header bytes  : %lu (%.2f%%)\n",
                st.bytes_posting_list_header,
                100.0 * double(st.bytes_posting_list_header) / double(total_bytes));
    std::printf("  posting list body bytes    : %lu (%.2f%%)\n",
                st.bytes_posting_list_body,
                100.0 * double(st.bytes_posting_list_body) / double(total_bytes));
}

bool walk_kix(const KixReader& kix, KixStats& st, uint64_t& total_bytes) {
    const uint8_t* base = kix.posting_file();
    const uint32_t tbl = kix.table_size();
    total_bytes = kix.posting_file_size();
    for (uint32_t kmer = 0; kmer < tbl; kmer++) {
        uint64_t off = kix.posting_list_offset(kmer);
        uint64_t len = kix.posting_list_byte_length(kmer);
        if (len == 0) continue;
        // Per-posting-list header is just [u32 distinct_count].
        if (len < 4) return false;
        uint32_t count;
        std::memcpy(&count, base + off, 4);
        const uint64_t body_bytes = len - 4;
        st.n_kmers_nonempty++;
        st.total_distinct_count += count;
        st.bytes_posting_list_header += 4;
        st.bytes_posting_list_body += body_bytes;
    }
    return true;
}

// Walk .kpx skipping k-mers whose .kix posting list is empty.  Empty
// k-mers in a production index have pos_offset == 0 (the builder does
// not emit a placeholder for them), which aliases the first non-empty
// k-mer's posting list — so we MUST consult the .kix to know which
// slots actually carry data.
bool walk_kpx(const KixReader& kix, const KpxReader& kpx,
              KpxStats& st, uint64_t& total_bytes) {
    const uint8_t* base = kpx.posting_file();
    total_bytes = kpx.posting_file_size();
    const uint32_t tbl = kpx.table_size();
    const uint8_t* end = base + total_bytes;
    uint64_t skipped_empty = 0;
    // walk_kpx_posting needs kix_count to size the kind map and derive
    // per-kind sub-bucket counts.  Read it from the leading u32 of
    // every .kix posting list.
    const uint8_t* kix_base = kix.posting_file();
    for (uint32_t kmer = 0; kmer < tbl; kmer++) {
        const uint64_t kix_len = kix.posting_list_byte_length(kmer);
        if (kix_len == 0) {
            skipped_empty++;
            continue;
        }
        const uint64_t kix_off = kix.posting_list_offset(kmer);
        uint32_t kix_count = 0;
        if (kix_len >= 4) {
            std::memcpy(&kix_count, kix_base + kix_off, 4);
        }
        uint64_t off = kpx.pos_offset(kmer);
        if (off > total_bytes) {
            std::fprintf(stderr, "kmer %u: offset %lu past end %lu\n",
                         kmer, off, total_bytes);
            return false;
        }
        if (!walk_kpx_posting(base + off, end, kix_count, st)) {
            std::fprintf(stderr, "kmer %u: walk failed at offset %lu\n", kmer, off);
            return false;
        }
    }
    std::printf("(skipped %lu empty kmers via .kix)\n", skipped_empty);
    // Re-derive total_position_count from the stream stats since
    // walk_kpx_posting accumulates per sub-bucket.
    st.total_position_count =
        st.partition_positions + st.short1_positions + st.short2_positions;
    return true;
}

void usage() {
    std::fprintf(stderr,
        "Usage: index_compress_stats --kix <path> [--kpx <path>]\n"
        "  .kix is always required.  When --kpx is given the .kix is used to\n"
        "  identify and skip empty-posting-list k-mers (which alias to offset 0\n"
        "  in the .kpx because the builder writes no placeholder).\n");
}

} // namespace

int main(int argc, char** argv) {
    std::string kix_path, kpx_path;
    for (int i = 1; i < argc; i++) {
        std::string s = argv[i];
        if (s == "--kix" && i + 1 < argc) kix_path = argv[++i];
        else if (s == "--kpx" && i + 1 < argc) kpx_path = argv[++i];
        else { usage(); return 1; }
    }
    if (kix_path.empty()) { usage(); return 1; }

    KixReader kix;
    if (!kix.open(kix_path)) {
        std::fprintf(stderr, "failed to open .kix: %s\n", kix_path.c_str());
        return 2;
    }
    std::printf("--- %s ---\n", kix_path.c_str());
    std::printf("k=%d t=%d num_seqs=%u total_distinct_postings=%lu\n",
                kix.k(), int(kix.t()), kix.num_sequences(),
                (unsigned long)kix.total_distinct_postings());
    {
        KixStats st;
        uint64_t total_bytes = 0;
        if (!walk_kix(kix, st, total_bytes)) {
            std::fprintf(stderr, "failed to walk .kix\n");
            return 3;
        }
        print_kix_report(st, total_bytes);
    }

    if (!kpx_path.empty()) {
        KpxReader kpx;
        if (!kpx.open(kpx_path)) {
            std::fprintf(stderr, "failed to open .kpx: %s\n", kpx_path.c_str());
            return 2;
        }
        if (kpx.k() != kix.k()) {
            std::fprintf(stderr, ".kix k=%d but .kpx k=%d mismatch\n",
                         kix.k(), kpx.k());
            return 3;
        }
        std::printf("\n--- %s ---\n", kpx_path.c_str());
        std::printf("k=%d t=%d total_position_count=%lu\n",
                    kpx.k(), int(kpx.t()),
                    (unsigned long)kpx.total_position_count());
        KpxStats st;
        uint64_t total_bytes = 0;
        if (!walk_kpx(kix, kpx, st, total_bytes)) {
            std::fprintf(stderr, "failed to walk .kpx\n");
            return 3;
        }
        print_kpx_report(st, total_bytes);
    }
    return 0;
}
