#include "index/index_builder.hpp"
#include "io/blastdb_reader.hpp"
#include "core/config.hpp"
#include "core/types.hpp"
#include "core/kmer_encoding.hpp"
#include "core/varint.hpp"
#include "index/ksx_writer.hpp"
#include "index/kix_format.hpp"
#include "index/kpx_format.hpp"
#include "util/logger.hpp"
#include "util/progress.hpp"

#include <cstdio>
#include <cstring>
#include <algorithm>
#include <numeric>
#include <vector>
#include <string>
#include <filesystem>

namespace ikafssn {

// Temporary entry for posting buffer.
struct TempEntry {
    uint32_t kmer_value;
    uint32_t seq_id;
    uint32_t pos;
};

static_assert(sizeof(TempEntry) == 12, "TempEntry must be 12 bytes");

// Determine partition for a kmer based on upper bits.
static inline uint32_t partition_of(uint32_t kmer, int partition_bits, int k) {
    if (partition_bits == 0) return 0;
    return (kmer >> (2 * k - partition_bits)) & ((1u << partition_bits) - 1);
}

// Compute ceiling log2 for powers of 2.
static inline int log2_ceil(int n) {
    int bits = 0;
    int v = n - 1;
    while (v > 0) { v >>= 1; bits++; }
    return bits;
}

template <typename KmerInt>
bool build_index(BlastDbReader& db,
                 const IndexBuilderConfig& config,
                 const std::string& output_prefix,
                 uint16_t volume_index,
                 uint16_t total_volumes,
                 const std::string& db_name,
                 const Logger& logger) {

    const int k = config.k;
    const uint64_t tbl_size = table_size(k);
    const uint32_t num_seqs = db.num_sequences();
    const int num_partitions = config.partitions;
    const int partition_bits = log2_ceil(num_partitions);

    logger.info("Building index: k=%d, sequences=%u, partitions=%d",
                k, num_seqs, num_partitions);

    // File paths (.tmp during construction)
    std::string ksx_tmp = output_prefix + ".ksx.tmp";
    std::string kix_tmp = output_prefix + ".kix.tmp";
    std::string kpx_tmp = output_prefix + ".kpx.tmp";
    std::string ksx_final = output_prefix + ".ksx";
    std::string kix_final = output_prefix + ".kix";
    std::string kpx_final = output_prefix + ".kpx";

    // =========== Phase 0: Metadata collection -> .ksx ===========
    logger.info("Phase 0: collecting metadata...");
    {
        KsxWriter ksx;
        Progress prog("Phase 0", num_seqs, config.verbose);
        for (uint32_t oid = 0; oid < num_seqs; oid++) {
            uint32_t slen = db.seq_length(oid);
            std::string acc = db.get_accession(oid);
            ksx.add_sequence(slen, acc);
            prog.update(oid + 1);
        }
        prog.finish();

        if (!ksx.write(ksx_tmp)) {
            logger.error("Failed to write %s", ksx_tmp.c_str());
            return false;
        }
        logger.info("Phase 0: wrote %s (%u sequences)", ksx_tmp.c_str(), num_seqs);
    }

    // =========== Phase 1: Counting pass ===========
    logger.info("Phase 1: counting k-mers...");
    std::vector<uint64_t> counts64(tbl_size, 0);
    {
        KmerScanner<KmerInt> scanner(k);
        Progress prog("Phase 1", num_seqs, config.verbose);

        for (uint32_t oid = 0; oid < num_seqs; oid++) {
            std::string seq = db.get_sequence(oid);
            scanner.scan(seq.data(), seq.size(),
                [&counts64](uint32_t /*pos*/, KmerInt kmer) {
                    counts64[kmer]++;
                });
            prog.update(oid + 1);
        }
        prog.finish();
    }

    // Convert uint64 -> uint32 with overflow check
    std::vector<uint32_t> counts(tbl_size);
    uint64_t total_postings = 0;
    for (uint64_t i = 0; i < tbl_size; i++) {
        if (counts64[i] > UINT32_MAX) {
            logger.error("k-mer %lu has count %lu which exceeds uint32_t. "
                         "Use a larger k value.", i, counts64[i]);
            std::remove(ksx_tmp.c_str());
            return false;
        }
        counts[i] = static_cast<uint32_t>(counts64[i]);
        total_postings += counts64[i];
    }
    counts64.clear();
    counts64.shrink_to_fit();

    logger.info("Phase 1: total postings = %lu", static_cast<unsigned long>(total_postings));

    // Apply max_freq_build exclusion
    if (config.max_freq_build > 0) {
        uint64_t excluded = 0;
        for (uint64_t i = 0; i < tbl_size; i++) {
            if (counts[i] > config.max_freq_build) {
                total_postings -= counts[i];
                counts[i] = 0;
                excluded++;
            }
        }
        logger.info("Phase 1: excluded %lu high-frequency k-mers (threshold=%lu)",
                    static_cast<unsigned long>(excluded),
                    static_cast<unsigned long>(config.max_freq_build));
    }

    // =========== Phase 2 & 3: Posting write with partition + buffer ===========
    // Method B: write postings first, record offsets, seek back to write tables.
    //
    // We open kix and kpx files, write headers and reserve table space,
    // then process partitions one by one.

    logger.info("Phase 2-3: writing postings...");

    // Open kix file
    FILE* kix_fp = std::fopen(kix_tmp.c_str(), "wb");
    if (!kix_fp) {
        logger.error("Cannot open %s for writing", kix_tmp.c_str());
        std::remove(ksx_tmp.c_str());
        return false;
    }

    // Open kpx file
    FILE* kpx_fp = std::fopen(kpx_tmp.c_str(), "wb");
    if (!kpx_fp) {
        logger.error("Cannot open %s for writing", kpx_tmp.c_str());
        std::fclose(kix_fp);
        std::remove(ksx_tmp.c_str());
        return false;
    }

    // Write kix header placeholder (will be overwritten at finalize)
    KixHeader kix_hdr{};
    std::fwrite(&kix_hdr, sizeof(kix_hdr), 1, kix_fp);

    // Reserve offsets table (uint64 * tbl_size)
    const uint64_t kix_offsets_pos = sizeof(KixHeader);
    std::vector<uint64_t> kix_offsets(tbl_size, 0);
    std::fwrite(kix_offsets.data(), sizeof(uint64_t), tbl_size, kix_fp);

    // Reserve counts table (uint32 * tbl_size) - we already have counts
    const uint64_t kix_counts_pos = kix_offsets_pos + sizeof(uint64_t) * tbl_size;
    std::fwrite(counts.data(), sizeof(uint32_t), tbl_size, kix_fp);

    // kix posting data starts here
    const uint64_t kix_posting_start = kix_counts_pos + sizeof(uint32_t) * tbl_size;

    // Write kpx header placeholder
    KpxHeader kpx_hdr{};
    std::fwrite(&kpx_hdr, sizeof(kpx_hdr), 1, kpx_fp);

    // Reserve pos_offsets table
    const uint64_t kpx_offsets_pos = sizeof(KpxHeader);
    std::vector<uint64_t> kpx_offsets(tbl_size, 0);
    std::fwrite(kpx_offsets.data(), sizeof(uint64_t), tbl_size, kpx_fp);

    // kpx posting data starts here
    const uint64_t kpx_posting_start = kpx_offsets_pos + sizeof(uint64_t) * tbl_size;

    // Current write positions in posting data (relative to posting start)
    uint64_t kix_data_pos = 0;
    uint64_t kpx_data_pos = 0;

    // Estimate buffer capacity
    uint64_t est_partition_postings = (total_postings + num_partitions - 1) / num_partitions;
    uint64_t buffer_entries_limit = config.buffer_size / sizeof(TempEntry);
    if (est_partition_postings > buffer_entries_limit) {
        logger.warn("Estimated partition size (%lu entries) exceeds buffer capacity (%lu entries). "
                    "Increase -buffer_size or -partitions.",
                    static_cast<unsigned long>(est_partition_postings),
                    static_cast<unsigned long>(buffer_entries_limit));
    }

    // Process each partition
    KmerScanner<KmerInt> scanner(k);
    uint8_t varint_buf[5];

    for (int p = 0; p < num_partitions; p++) {
        logger.info("  Partition %d/%d...", p + 1, num_partitions);

        // Collect all entries for this partition into buffer
        std::vector<TempEntry> buffer;
        buffer.reserve(std::min(est_partition_postings, buffer_entries_limit));

        Progress prog("  Partition scan", num_seqs, config.verbose);
        for (uint32_t oid = 0; oid < num_seqs; oid++) {
            std::string seq = db.get_sequence(oid);
            scanner.scan(seq.data(), seq.size(),
                [&](uint32_t pos, KmerInt kmer) {
                    uint32_t kval = static_cast<uint32_t>(kmer);
                    if (counts[kval] == 0) return; // excluded
                    if (static_cast<int>(partition_of(kval, partition_bits, k)) != p) return;
                    buffer.push_back({kval, oid, pos});
                });
            prog.update(oid + 1);
        }
        prog.finish();

        if (buffer.empty()) continue;

        // Sort by kmer, then seq_id, then pos
        std::sort(buffer.begin(), buffer.end(),
            [](const TempEntry& a, const TempEntry& b) {
                if (a.kmer_value != b.kmer_value) return a.kmer_value < b.kmer_value;
                if (a.seq_id != b.seq_id) return a.seq_id < b.seq_id;
                return a.pos < b.pos;
            });

        // Write sorted postings grouped by kmer
        size_t i = 0;
        while (i < buffer.size()) {
            uint32_t cur_kmer = buffer[i].kmer_value;

            // Find range [i, j) for this kmer
            size_t j = i;
            while (j < buffer.size() && buffer[j].kmer_value == cur_kmer) {
                j++;
            }

            // Record offsets for this kmer
            kix_offsets[cur_kmer] = kix_data_pos;
            kpx_offsets[cur_kmer] = kpx_data_pos;

            // Write delta-compressed ID postings to kix
            {
                uint32_t prev_id = 0;
                for (size_t e = i; e < j; e++) {
                    uint32_t delta = (e == i) ? buffer[e].seq_id
                                              : buffer[e].seq_id - prev_id;
                    prev_id = buffer[e].seq_id;
                    size_t n = varint_encode(delta, varint_buf);
                    std::fwrite(varint_buf, 1, n, kix_fp);
                    kix_data_pos += n;
                }
            }

            // Write delta-compressed pos postings to kpx
            // Position deltas reset at sequence boundaries
            {
                uint32_t prev_id = UINT32_MAX; // force "new seq" on first
                uint32_t prev_pos = 0;
                for (size_t e = i; e < j; e++) {
                    bool new_seq = (buffer[e].seq_id != prev_id);
                    uint32_t val;
                    if (new_seq) {
                        val = buffer[e].pos; // raw position (delta reset)
                    } else {
                        val = buffer[e].pos - prev_pos; // delta
                    }
                    prev_id = buffer[e].seq_id;
                    prev_pos = buffer[e].pos;
                    size_t n = varint_encode(val, varint_buf);
                    std::fwrite(varint_buf, 1, n, kpx_fp);
                    kpx_data_pos += n;
                }
            }

            i = j;
        }

        logger.debug("  Partition %d: %lu entries written", p + 1,
                     static_cast<unsigned long>(buffer.size()));
    }

    // =========== Phase 4: Finalize ===========
    logger.info("Phase 4: finalizing...");

    // Write kix header
    std::memcpy(kix_hdr.magic, KIX_MAGIC, 4);
    kix_hdr.format_version = KIX_FORMAT_VERSION;
    kix_hdr.k = static_cast<uint8_t>(k);
    kix_hdr.kmer_type = kmer_type_for_k(k);
    kix_hdr.num_sequences = num_seqs;
    kix_hdr.total_postings = total_postings;
    kix_hdr.flags = KIX_FLAG_HAS_KSX;
    kix_hdr.volume_index = volume_index;
    kix_hdr.total_volumes = total_volumes;
    size_t name_len = std::min(db_name.size(), size_t(32));
    kix_hdr.db_name_len = static_cast<uint16_t>(name_len);
    std::memcpy(kix_hdr.db_name, db_name.c_str(), name_len);

    std::fseek(kix_fp, 0, SEEK_SET);
    std::fwrite(&kix_hdr, sizeof(kix_hdr), 1, kix_fp);

    // Write kix offsets table
    std::fwrite(kix_offsets.data(), sizeof(uint64_t), tbl_size, kix_fp);

    // Counts table was already written correctly during reservation
    // (we wrote the actual counts array, not zeros)

    std::fclose(kix_fp);

    // Write kpx header
    std::memcpy(kpx_hdr.magic, KPX_MAGIC, 4);
    kpx_hdr.format_version = KPX_FORMAT_VERSION;
    kpx_hdr.k = static_cast<uint8_t>(k);
    kpx_hdr.total_postings = total_postings;

    std::fseek(kpx_fp, 0, SEEK_SET);
    std::fwrite(&kpx_hdr, sizeof(kpx_hdr), 1, kpx_fp);

    // Write kpx offsets table
    std::fwrite(kpx_offsets.data(), sizeof(uint64_t), tbl_size, kpx_fp);

    std::fclose(kpx_fp);

    // Rename .tmp files to final names
    if (std::rename(ksx_tmp.c_str(), ksx_final.c_str()) != 0) {
        logger.error("Failed to rename %s -> %s", ksx_tmp.c_str(), ksx_final.c_str());
        return false;
    }
    if (std::rename(kix_tmp.c_str(), kix_final.c_str()) != 0) {
        logger.error("Failed to rename %s -> %s", kix_tmp.c_str(), kix_final.c_str());
        return false;
    }
    if (std::rename(kpx_tmp.c_str(), kpx_final.c_str()) != 0) {
        logger.error("Failed to rename %s -> %s", kpx_tmp.c_str(), kpx_final.c_str());
        return false;
    }

    logger.info("Index built: %s (.kix, .kpx, .ksx)", output_prefix.c_str());
    return true;
}

// Explicit template instantiations
template bool build_index<uint16_t>(BlastDbReader&, const IndexBuilderConfig&,
    const std::string&, uint16_t, uint16_t, const std::string&, const Logger&);
template bool build_index<uint32_t>(BlastDbReader&, const IndexBuilderConfig&,
    const std::string&, uint16_t, uint16_t, const std::string&, const Logger&);

} // namespace ikafssn
