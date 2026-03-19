#include "index/index_builder.hpp"
#include "io/blastdb_reader.hpp"
#include "core/config.hpp"
#include "core/types.hpp"
#include "core/kmer_encoding.hpp"
#include "core/packed_kmer_scanner.hpp"
#include "core/ambiguity_parser.hpp"
#include "core/varint.hpp"
#include "core/spaced_seed.hpp"
#include "index/ksx_writer.hpp"
#include "index/kix_format.hpp"
#include "index/kpx_format.hpp"
#include "util/logger.hpp"
#include "util/progress.hpp"

#include <cstdio>
#include <cstring>
#include <algorithm>
#include <chrono>
#include <numeric>
#include <vector>
#include <string>
#include <filesystem>

#include <tbb/parallel_for.h>
#include <tbb/parallel_sort.h>
#include <tbb/blocked_range.h>
#include <tbb/combinable.h>

#include <atomic>

namespace ikafssn {

// Temporary entry for posting buffer.
struct TempEntry {
    uint32_t kmer_value;
    uint32_t seq_id;
    uint32_t pos;
};

static_assert(sizeof(TempEntry) == 12, "TempEntry must be 12 bytes");

// Determine partition for a kmer based on upper bits.
// effective_bits is the total number of significant bits (2*k, or 2*k+1 for tagged).
static inline uint32_t partition_of(uint32_t kmer, int partition_bits, int effective_bits) {
    if (partition_bits == 0) return 0;
    return (kmer >> (effective_bits - partition_bits)) & ((1u << partition_bits) - 1);
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
    const uint32_t num_seqs = db.num_sequences();

    logger.info("Building index: k=%d, sequences=%u", k, num_seqs);

    // File paths (.tmp during construction, renamed to final on success)
    std::string ksx_tmp = output_prefix + ".ksx.tmp";
    std::string kix_tmp = output_prefix + ".kix.tmp";
    std::string ksx_final = output_prefix + ".ksx";
    std::string kix_final = output_prefix + ".kix";
    std::string kpx_tmp, kpx_final;
    if (!config.skip_kpx) {
        kpx_tmp = output_prefix + ".kpx.tmp";
        kpx_final = output_prefix + ".kpx";
    }

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

    // Pre-compute spaced seed masks and tags (shared across all phases).
    // For "both" mode, tags embed a 1-bit template identity into k-mer values.
    std::vector<uint32_t> seed_masks;
    std::vector<KmerInt> seed_tags; // pre-shifted tag values; empty = no tagging
    if (config.t > 0) {
        auto [m, t] = get_tagged_masks<KmerInt>(k, config.t,
                           config.template_type, config.template_type);
        seed_masks = std::move(m);
        seed_tags = std::move(t);
    }

    // Effective table size: doubles for "both" mode (tag bit per mask)
    const int num_masks = static_cast<int>(seed_masks.size());
    const uint32_t tbl_size = spaced_table_size(k, num_masks > 0 ? num_masks : 1);
    // Effective bit width for partitioning (2*k, or 2*k+1 for tagged both mode)
    const int effective_bits = 2 * k + ((num_masks > 1) ? 1 : 0);

    // =========== Phase 1: Counting pass (TBB parallel) ===========
    logger.info("Phase 1: counting k-mers (threads=%d)...", config.threads);
    std::vector<uint64_t> counts64(tbl_size, 0);
    {
        // Always use parallel path; TBB respects global_control parallelism
        // (threads==1 degrades gracefully to sequential execution)
        tbb::combinable<std::vector<uint64_t>> local_counts(
            [&tbl_size]() { return std::vector<uint64_t>(tbl_size, 0); });

        tbb::parallel_for(
            tbb::blocked_range<uint32_t>(0, num_seqs, 64),
            [&](const tbb::blocked_range<uint32_t>& range) {
                auto& my_counts = local_counts.local();
                PackedKmerScanner<KmerInt> scanner(k);
                for (uint32_t oid = range.begin(); oid < range.end(); oid++) {
                    auto raw = db.get_raw_sequence(oid);
                    auto ambig = AmbiguityParser::parse(raw.ambig_data, raw.ambig_bytes);
                    if (config.t > 0) {
                        scanner.scan_spaced(raw.ncbi2na_data, raw.seq_length, ambig,
                            seed_masks, static_cast<int>(config.t),
                            [&my_counts](uint32_t /*pos*/, KmerInt kmer) {
                                my_counts[kmer]++;
                            },
                            [&my_counts](uint32_t /*pos*/, KmerInt base_kmer,
                                         const AmbigInfo* infos, int count) {
                                expand_ambig_kmer_multi<KmerInt>(base_kmer, infos, count,
                                    [&my_counts](KmerInt expanded) {
                                        my_counts[expanded]++;
                                    });
                            },
                            config.max_degen_expand,
                            seed_tags);
                    } else {
                        scanner.scan(raw.ncbi2na_data, raw.seq_length, ambig,
                            [&my_counts](uint32_t /*pos*/, KmerInt kmer) {
                                my_counts[kmer]++;
                            },
                            [&my_counts](uint32_t /*pos*/, KmerInt base_kmer,
                                         const AmbigInfo* infos, int count) {
                                expand_ambig_kmer_multi<KmerInt>(base_kmer, infos, count,
                                    [&my_counts](KmerInt expanded) {
                                        my_counts[expanded]++;
                                    });
                            },
                            config.max_degen_expand);
                    }
                    db.ret_raw_sequence(raw);
                }
            });

        // Reduce thread-local counts
        local_counts.combine_each([&counts64, &tbl_size](const std::vector<uint64_t>& lc) {
            for (uint32_t i = 0; i < tbl_size; i++) {
                counts64[i] += lc[i];
            }
        });
    }

    // Convert uint64 -> uint32 with overflow check
    std::vector<uint32_t> counts(tbl_size);
    uint64_t total_postings = 0;
    for (uint32_t i = 0; i < tbl_size; i++) {
        if (counts64[i] > UINT32_MAX) {
            logger.error("k-mer %u has count %lu which exceeds uint32_t. "
                         "Use a larger k value.", i, counts64[i]);
            std::remove(ksx_tmp.c_str());
            return false;
        }
        counts[i] = static_cast<uint32_t>(counts64[i]);
        total_postings += counts64[i];
    }
    counts64.clear();
    counts64.shrink_to_fit();

    logger.info("Phase 1: total postings = %lu (estimated)", static_cast<unsigned long>(total_postings));

    // When using spaced seeds with both templates, two different masks may
    // produce the same (kmer, seq_id, pos) entry. Track actual counts after
    // deduplication in Phase 2-3 and rewrite the counts table in Phase 4.
    const bool need_dedup = (config.t > 0);
    std::vector<uint32_t> actual_counts(tbl_size, 0);
    uint64_t actual_total_postings = 0;

    // =========== Determine partition count from memory_limit ===========
    int num_partitions = 1;
    if (total_postings > 0) {
        uint64_t entries_limit = config.memory_limit / sizeof(TempEntry);
        while (static_cast<uint64_t>((total_postings + num_partitions - 1) / num_partitions)
               > entries_limit) {
            num_partitions *= 2;
        }
    }
    const int partition_bits = log2_ceil(num_partitions);

    if (config.memory_limit >= (uint64_t(1) << 30))
        logger.info("Phase 2-3: writing postings (partitions=%d, memory_limit=%luG)...",
                    num_partitions,
                    static_cast<unsigned long>(config.memory_limit >> 30));
    else
        logger.info("Phase 2-3: writing postings (partitions=%d, memory_limit=%luM)...",
                    num_partitions,
                    static_cast<unsigned long>(config.memory_limit >> 20));

    // Open kix file
    FILE* kix_fp = std::fopen(kix_tmp.c_str(), "wb");
    if (!kix_fp) {
        logger.error("Cannot open %s for writing", kix_tmp.c_str());
        std::remove(ksx_tmp.c_str());
        return false;
    }

    // Open kpx file (skip if mode 1)
    FILE* kpx_fp = nullptr;
    if (!config.skip_kpx) {
        kpx_fp = std::fopen(kpx_tmp.c_str(), "wb");
        if (!kpx_fp) {
            logger.error("Cannot open %s for writing", kpx_tmp.c_str());
            std::fclose(kix_fp);
            std::remove(ksx_tmp.c_str());
            return false;
        }
    }

    // Write kix header placeholder (will be overwritten at finalize)
    KixHeader kix_hdr{};
    std::fwrite(&kix_hdr, sizeof(kix_hdr), 1, kix_fp);

    // Reserve offsets table (always uint64 in builder; filter may compact to uint32)
    const uint64_t kix_offsets_pos = sizeof(KixHeader);
    std::vector<uint64_t> kix_offsets(tbl_size + 1, 0);
    std::fwrite(kix_offsets.data(), sizeof(uint64_t), tbl_size + 1, kix_fp);

    // kix posting data starts here (no counts table in format v3)
    const uint64_t kix_posting_start = kix_offsets_pos + sizeof(uint64_t) * (tbl_size + 1);

    // Write kpx header placeholder (skip if mode 1)
    KpxHeader kpx_hdr{};
    std::vector<uint64_t> kpx_offsets;
    if (!config.skip_kpx) {
        std::fwrite(&kpx_hdr, sizeof(kpx_hdr), 1, kpx_fp);

        // Reserve pos_offsets table
        kpx_offsets.resize(tbl_size, 0);
        std::fwrite(kpx_offsets.data(), sizeof(uint64_t), tbl_size, kpx_fp);
    }

    // Current write positions in posting data (relative to posting start)
    uint64_t kix_data_pos = 0;
    uint64_t kpx_data_pos = 0;

    uint64_t est_partition_postings = (total_postings + num_partitions - 1) / num_partitions;
    uint64_t reserve_entries = config.memory_limit / sizeof(TempEntry);

    // Process each partition (sequentially to respect memory constraints)
    uint8_t varint_buf[5];

    for (int p = 0; p < num_partitions; p++) {
        logger.info("  Partition %d/%d...", p + 1, num_partitions);

        // Parallel scan: collect entries for this partition using thread-local buffers
        tbb::combinable<std::vector<TempEntry>> local_buffers;

        std::atomic<uint32_t> progress_counter{0};
        auto progress_start = std::chrono::steady_clock::now();

        tbb::parallel_for(
            tbb::blocked_range<uint32_t>(0, num_seqs, 64),
            [&](const tbb::blocked_range<uint32_t>& range) {
                auto& my_buffer = local_buffers.local();
                PackedKmerScanner<KmerInt> scanner(k);
                for (uint32_t oid = range.begin(); oid < range.end(); oid++) {
                    auto raw = db.get_raw_sequence(oid);
                    auto ambig = AmbiguityParser::parse(raw.ambig_data, raw.ambig_bytes);
                    auto normal_cb = [&](uint32_t pos, KmerInt kmer) {
                        uint32_t kval = static_cast<uint32_t>(kmer);
                        if (counts[kval] == 0) return;
                        if (static_cast<int>(partition_of(kval, partition_bits, effective_bits)) != p) return;
                        my_buffer.push_back({kval, oid, pos});
                    };
                    auto ambig_cb = [&](uint32_t pos, KmerInt base_kmer,
                                        const AmbigInfo* infos, int count) {
                        expand_ambig_kmer_multi<KmerInt>(base_kmer, infos, count,
                            [&](KmerInt expanded) {
                                uint32_t kval = static_cast<uint32_t>(expanded);
                                if (counts[kval] == 0) return;
                                if (static_cast<int>(partition_of(kval, partition_bits, effective_bits)) != p) return;
                                my_buffer.push_back({kval, oid, pos});
                            });
                    };
                    if (config.t > 0) {
                        scanner.scan_spaced(raw.ncbi2na_data, raw.seq_length, ambig,
                            seed_masks, static_cast<int>(config.t),
                            normal_cb, ambig_cb,
                            config.max_degen_expand,
                            seed_tags);
                    } else {
                        scanner.scan(raw.ncbi2na_data, raw.seq_length, ambig,
                            normal_cb, ambig_cb,
                            config.max_degen_expand);
                    }
                    db.ret_raw_sequence(raw);
                }
                // Atomic progress update (coarse-grained per chunk)
                uint32_t done = progress_counter.fetch_add(
                    range.end() - range.begin(), std::memory_order_relaxed)
                    + (range.end() - range.begin());
                if (config.verbose && done % 1000 < (range.end() - range.begin())) {
                    auto now = std::chrono::steady_clock::now();
                    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(
                        now - progress_start).count();
                    std::fprintf(stderr, "\r  Partition scan: %.1f%% (%u/%u) [%lds]",
                                 100.0 * done / num_seqs, done, num_seqs,
                                 static_cast<long>(elapsed));
                    std::fflush(stderr);
                }
            });

        if (config.verbose) {
            auto now = std::chrono::steady_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(
                now - progress_start).count();
            std::fprintf(stderr, "\r  Partition scan: done (%u items, %lds)\n",
                         num_seqs, static_cast<long>(elapsed));
            std::fflush(stderr);
        }

        // Merge thread-local buffers into single buffer
        std::vector<TempEntry> buffer;
        buffer.reserve(std::min(est_partition_postings, reserve_entries));
        local_buffers.combine_each([&buffer](std::vector<TempEntry>& local) {
            buffer.insert(buffer.end(),
                          std::make_move_iterator(local.begin()),
                          std::make_move_iterator(local.end()));
            // Release local memory
            std::vector<TempEntry>().swap(local);
        });

        if (buffer.empty()) continue;

        // Parallel sort by kmer, then seq_id, then pos
        tbb::parallel_sort(buffer.begin(), buffer.end(),
            [](const TempEntry& a, const TempEntry& b) {
                if (a.kmer_value != b.kmer_value) return a.kmer_value < b.kmer_value;
                if (a.seq_id != b.seq_id) return a.seq_id < b.seq_id;
                return a.pos < b.pos;
            });

        // Deduplicate entries with identical (kmer_value, seq_id, pos).
        // This occurs when multiple spaced seed masks produce the same k-mer
        // at the same position for the same sequence.
        if (need_dedup) {
            auto new_end = std::unique(buffer.begin(), buffer.end(),
                [](const TempEntry& a, const TempEntry& b) {
                    return a.kmer_value == b.kmer_value
                        && a.seq_id == b.seq_id
                        && a.pos == b.pos;
                });
            buffer.erase(new_end, buffer.end());
        }

        // Write sorted postings grouped by kmer (sequential — I/O bound)
        size_t i = 0;
        while (i < buffer.size()) {
            uint32_t cur_kmer = buffer[i].kmer_value;

            // Find range [i, j) for this kmer
            size_t j = i;
            while (j < buffer.size() && buffer[j].kmer_value == cur_kmer) {
                j++;
            }

            // Record offsets and actual counts for this kmer
            kix_offsets[cur_kmer] = kix_data_pos;
            if (!config.skip_kpx) kpx_offsets[cur_kmer] = kpx_data_pos;
            actual_counts[cur_kmer] = static_cast<uint32_t>(j - i);
            actual_total_postings += (j - i);

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

            // Write delta-compressed pos postings to kpx (skip if mode 1)
            if (!config.skip_kpx) {
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

    if (need_dedup && actual_total_postings < total_postings) {
        logger.info("Deduplicated: %lu -> %lu postings (removed %lu duplicates)",
                    static_cast<unsigned long>(total_postings),
                    static_cast<unsigned long>(actual_total_postings),
                    static_cast<unsigned long>(total_postings - actual_total_postings));
    }

    // Forward-fill kix_offsets: empty k-mers get the same offset as the next
    // non-empty k-mer (or the sentinel). This ensures offsets[i+1]-offsets[i]==0
    // for empty k-mers.
    {
        uint64_t fill = kix_data_pos; // sentinel value for trailing empties
        for (int32_t i = static_cast<int32_t>(tbl_size) - 1; i >= 0; i--) {
            if (actual_counts[i] > 0) {
                fill = kix_offsets[i];
            } else {
                kix_offsets[i] = fill;
            }
        }
    }

    // =========== Phase 4: Finalize ===========
    logger.info("Phase 4: finalizing...");

    // Set sentinel offset
    kix_offsets[tbl_size] = kix_data_pos;

    // Determine offset types based on posting data sizes
    const bool kix_offset32 = (kix_data_pos <= UINT32_MAX);
    const bool kpx_offset32 = (!config.skip_kpx && kpx_data_pos <= UINT32_MAX);

    // Close the in-progress kix file; we'll rewrite it with the correct layout.
    std::fclose(kix_fp);

    // Rewrite .kix with correct offset width
    {
        // Read the posting data from the temp file (skip header + uint64 offsets)
        FILE* rd = std::fopen(kix_tmp.c_str(), "rb");
        std::fseek(rd, static_cast<long>(kix_posting_start), SEEK_SET);
        std::vector<uint8_t> posting_blob(kix_data_pos);
        if (kix_data_pos > 0) {
            std::fread(posting_blob.data(), 1, kix_data_pos, rd);
        }
        std::fclose(rd);

        // Write the final kix file
        FILE* wr = std::fopen(kix_tmp.c_str(), "wb");

        std::memcpy(kix_hdr.magic, KIX_MAGIC, 4);
        kix_hdr.format_version = KIX_FORMAT_VERSION;
        kix_hdr.k = static_cast<uint8_t>(k);
        kix_hdr.kmer_type = kmer_type_for(k, config.t);
        kix_hdr.num_sequences = num_seqs;
        kix_hdr.total_postings = actual_total_postings;
        kix_hdr.flags = KIX_FLAG_HAS_KSX | (kix_offset32 ? KIX_FLAG_OFFSET32 : 0);
        kix_hdr.volume_index = volume_index;
        kix_hdr.total_volumes = total_volumes;
        size_t name_len = std::min(db_name.size(), size_t(32));
        kix_hdr.db_len = static_cast<uint16_t>(name_len);
        std::memcpy(kix_hdr.db, db_name.c_str(), name_len);
        kix_hdr.t = config.t;
        kix_hdr.template_type = config.template_type;
        std::fwrite(&kix_hdr, sizeof(kix_hdr), 1, wr);

        if (kix_offset32) {
            std::vector<uint32_t> off32(tbl_size + 1);
            for (uint32_t i = 0; i <= tbl_size; i++)
                off32[i] = static_cast<uint32_t>(kix_offsets[i]);
            std::fwrite(off32.data(), sizeof(uint32_t), tbl_size + 1, wr);
        } else {
            std::fwrite(kix_offsets.data(), sizeof(uint64_t), tbl_size + 1, wr);
        }

        if (!posting_blob.empty()) {
            std::fwrite(posting_blob.data(), 1, posting_blob.size(), wr);
        }
        std::fclose(wr);
    }

    // Rewrite .kpx with correct offset width (skip if mode 1)
    if (!config.skip_kpx) {
        // kpx posting start: header + uint64 offsets
        const uint64_t kpx_posting_start_pos = sizeof(KpxHeader) + sizeof(uint64_t) * tbl_size;

        FILE* rd = std::fopen(kpx_tmp.c_str(), "rb");
        std::fseek(rd, static_cast<long>(kpx_posting_start_pos), SEEK_SET);
        std::vector<uint8_t> posting_blob(kpx_data_pos);
        if (kpx_data_pos > 0) {
            std::fread(posting_blob.data(), 1, kpx_data_pos, rd);
        }
        std::fclose(rd);

        FILE* wr = std::fopen(kpx_tmp.c_str(), "wb");

        std::memcpy(kpx_hdr.magic, KPX_MAGIC, 4);
        kpx_hdr.format_version = KPX_FORMAT_VERSION;
        kpx_hdr.k = static_cast<uint8_t>(k);
        kpx_hdr.t = config.t;
        kpx_hdr.template_type = config.template_type;
        kpx_hdr.total_postings = actual_total_postings;
        kpx_hdr.offset_type = kpx_offset32 ? 0 : 1;
        std::fwrite(&kpx_hdr, sizeof(kpx_hdr), 1, wr);

        if (kpx_offset32) {
            std::vector<uint32_t> off32(tbl_size);
            for (uint32_t i = 0; i < tbl_size; i++)
                off32[i] = static_cast<uint32_t>(kpx_offsets[i]);
            std::fwrite(off32.data(), sizeof(uint32_t), tbl_size, wr);
        } else {
            std::fwrite(kpx_offsets.data(), sizeof(uint64_t), tbl_size, wr);
        }

        if (!posting_blob.empty()) {
            std::fwrite(posting_blob.data(), 1, posting_blob.size(), wr);
        }
        std::fclose(wr);
    }

    // Rename .tmp files to final names (unless keep_tmp is set)
    if (!config.keep_tmp) {
        if (std::rename(ksx_tmp.c_str(), ksx_final.c_str()) != 0) {
            logger.error("Failed to rename %s -> %s", ksx_tmp.c_str(), ksx_final.c_str());
            return false;
        }
        if (std::rename(kix_tmp.c_str(), kix_final.c_str()) != 0) {
            logger.error("Failed to rename %s -> %s", kix_tmp.c_str(), kix_final.c_str());
            return false;
        }
        if (!config.skip_kpx) {
            if (std::rename(kpx_tmp.c_str(), kpx_final.c_str()) != 0) {
                logger.error("Failed to rename %s -> %s", kpx_tmp.c_str(), kpx_final.c_str());
                return false;
            }
        }
    }

    logger.info("Index built: %s (.kix%s, .ksx%s)", output_prefix.c_str(),
                config.skip_kpx ? "" : ", .kpx",
                config.keep_tmp ? " [tmp]" : "");
    return true;
}

// Explicit template instantiations
template bool build_index<uint16_t>(BlastDbReader&, const IndexBuilderConfig&,
    const std::string&, uint16_t, uint16_t, const std::string&, const Logger&);
template bool build_index<uint32_t>(BlastDbReader&, const IndexBuilderConfig&,
    const std::string&, uint16_t, uint16_t, const std::string&, const Logger&);

} // namespace ikafssn
