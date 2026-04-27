// Phase 5a — dump v3 .kix/.kpx postings to TSV for compression analysis.
//
// Throwaway tool: reads existing v3 (.kix v3 sentinel-offset + LEB128 deltas,
// .kpx v3 with within-seq position deltas) and emits TSV lines:
//
//     kmer_value \t count \t seq_id_csv \t pos_csv
//
// with absolute seq_id (no delta) and absolute position (no delta), one row
// per non-empty k-mer. The output feeds tools/analyze_postings.py which
// estimates v4 (FastPFOR* + Simple-8b) compressed size.
//
// Usage:
//     dump_postings --kix <path> [--kpx <path>] [--max-kmers N] [--min-count N]
//                   [--out <path>]
//
// If --kpx is omitted, only the .kix portion is dumped (pos column = empty).
// --max-kmers / --min-count gate the output to keep TSVs manageable.

#include "core/varint.hpp"
#include "core/config.hpp"
#include "index/kix_reader.hpp"
#include "index/kpx_reader.hpp"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

using namespace ikafssn;

namespace {

struct Args {
    std::string kix_path;
    std::string kpx_path;
    std::string out_path;
    uint64_t    max_kmers = 0;     // 0 = unlimited
    uint32_t    min_count = 1;     // skip empty k-mers by default
    uint32_t    stride = 1;        // emit every Nth non-empty k-mer (uniform sampling)
};

void usage() {
    std::fprintf(stderr,
        "Usage: dump_postings --kix <path> [--kpx <path>] [options]\n"
        "  --kix       <path>   v3 .kix file (required)\n"
        "  --kpx       <path>   v3 .kpx file (optional)\n"
        "  --out       <path>   output TSV (default stdout)\n"
        "  --max-kmers N        emit at most N non-empty k-mers (0 = all)\n"
        "  --min-count N        skip k-mers with count < N (default 1)\n"
        "  --stride N           sample every Nth non-empty k-mer (default 1)\n");
}

bool parse_args(int argc, char** argv, Args& a) {
    for (int i = 1; i < argc; i++) {
        std::string s = argv[i];
        auto need = [&](const char* name) -> const char* {
            if (i + 1 >= argc) {
                std::fprintf(stderr, "missing value for %s\n", name);
                return nullptr;
            }
            return argv[++i];
        };
        if (s == "--kix") { auto v = need("--kix"); if (!v) return false; a.kix_path = v; }
        else if (s == "--kpx") { auto v = need("--kpx"); if (!v) return false; a.kpx_path = v; }
        else if (s == "--out") { auto v = need("--out"); if (!v) return false; a.out_path = v; }
        else if (s == "--max-kmers") { auto v = need("--max-kmers"); if (!v) return false; a.max_kmers = std::strtoull(v, nullptr, 10); }
        else if (s == "--min-count") { auto v = need("--min-count"); if (!v) return false; a.min_count = static_cast<uint32_t>(std::strtoul(v, nullptr, 10)); }
        else if (s == "--stride") { auto v = need("--stride"); if (!v) return false; a.stride = static_cast<uint32_t>(std::strtoul(v, nullptr, 10)); if (a.stride == 0) a.stride = 1; }
        else if (s == "-h" || s == "--help") { usage(); return false; }
        else { std::fprintf(stderr, "unknown arg: %s\n", s.c_str()); usage(); return false; }
    }
    if (a.kix_path.empty()) { usage(); return false; }
    return true;
}

} // namespace

int main(int argc, char** argv) {
    Args args;
    if (!parse_args(argc, argv, args)) return 1;

    KixReader kix;
    if (!kix.open(args.kix_path)) {
        std::fprintf(stderr, "failed to open .kix: %s\n", args.kix_path.c_str());
        return 2;
    }

    KpxReader kpx;
    bool have_kpx = !args.kpx_path.empty();
    if (have_kpx && !kpx.open(args.kpx_path)) {
        std::fprintf(stderr, "failed to open .kpx: %s\n", args.kpx_path.c_str());
        return 2;
    }

    if (have_kpx) {
        if (kpx.k() != kix.k()) {
            std::fprintf(stderr, ".kix k=%d but .kpx k=%d mismatch\n",
                         kix.k(), kpx.k());
            return 3;
        }
    }

    FILE* out = stdout;
    if (!args.out_path.empty()) {
        out = std::fopen(args.out_path.c_str(), "w");
        if (!out) {
            std::fprintf(stderr, "failed to open output: %s\n", args.out_path.c_str());
            return 2;
        }
    }

    // TSV header (commented so awk users can skip with grep -v '^#')
    std::fprintf(out, "# k=%d t=%d num_seqs=%u total_postings=%lu\n",
                 kix.k(), (int)kix.t(), kix.num_sequences(),
                 (unsigned long)kix.total_postings());
    std::fprintf(out, "# kmer\tcount\tseq_ids\tpositions\n");

    const uint8_t* kix_data = kix.posting_data();
    const uint8_t* kpx_data = have_kpx ? kpx.posting_data() : nullptr;

    const uint32_t tbl = kix.table_size();
    uint64_t emitted = 0;
    uint64_t non_empty_seen = 0;

    std::vector<uint32_t> seq_ids;
    std::vector<uint32_t> positions;
    seq_ids.reserve(1024);
    positions.reserve(1024);

    for (uint32_t kmer = 0; kmer < tbl; kmer++) {
        if (args.max_kmers > 0 && emitted >= args.max_kmers) break;

        uint64_t kix_off  = kix.posting_offset(kmer);
        uint64_t kix_len  = kix.posting_byte_length(kmer);
        if (kix_len == 0) continue;

        // Stride-based uniform sampling over non-empty k-mers.
        if ((non_empty_seen++ % args.stride) != 0) continue;

        // Decode seq_ids (delta).
        seq_ids.clear();
        {
            const uint8_t* p = kix_data + kix_off;
            const uint8_t* e = p + kix_len;
            uint32_t prev = 0;
            bool first = true;
            while (p < e) {
                uint32_t d;
                p += varint_decode(p, d);
                uint32_t sid = first ? d : prev + d;
                seq_ids.push_back(sid);
                prev = sid;
                first = false;
            }
        }
        uint32_t cnt = static_cast<uint32_t>(seq_ids.size());
        if (cnt < args.min_count) continue;

        // Decode positions (delta within seq, reset on new seq).
        positions.clear();
        if (have_kpx) {
            uint64_t kpx_off_lo = kpx.pos_offset(kmer);
            // .kpx v3 has no sentinel; derive next offset by scanning forward.
            // We already know the count, decode exactly cnt varints.
            const uint8_t* p = kpx_data + kpx_off_lo;
            // Compute the delta-reset bitmask from seq_id deltas (first is reset).
            uint32_t prev_sid = 0;
            uint32_t prev_pos = 0;
            for (uint32_t i = 0; i < cnt; i++) {
                uint32_t v;
                p += varint_decode(p, v);
                bool new_seq = (i == 0) || (seq_ids[i] != prev_sid);
                uint32_t abs_pos = new_seq ? v : prev_pos + v;
                positions.push_back(abs_pos);
                prev_sid = seq_ids[i];
                prev_pos = abs_pos;
            }
        }

        // Emit row.
        std::fprintf(out, "%u\t%u\t", kmer, cnt);
        for (uint32_t i = 0; i < cnt; i++) {
            std::fprintf(out, "%u%s", seq_ids[i], (i + 1 == cnt ? "\t" : ","));
        }
        if (have_kpx) {
            for (uint32_t i = 0; i < cnt; i++) {
                std::fprintf(out, "%u%s", positions[i], (i + 1 == cnt ? "" : ","));
            }
        }
        std::fputc('\n', out);
        emitted++;
    }

    if (out != stdout) std::fclose(out);
    std::fprintf(stderr, "dump_postings: emitted %lu k-mers\n",
                 (unsigned long)emitted);
    return 0;
}
