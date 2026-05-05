#include "index/index_builder.hpp"
#include "io/blastdb_reader.hpp"
#include "core/config.hpp"
#include "core/types.hpp"
#include "core/kmer_encoding.hpp"
#include "core/packed_kmer_scanner.hpp"
#include "core/ambiguity_parser.hpp"
#include "core/spaced_seed.hpp"
#include "index/ksx_writer.hpp"
#include "index/kix_format.hpp"
#include "index/dictionary_io.hpp"
#include "index/kpx_format.hpp"
#include "index/pfd_codec.hpp"
#include "index/seq_id_dedup.hpp"
#include "index/parallel_sort_dispatch.hpp"
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
#include <tbb/blocked_range.h>
#include <tbb/combinable.h>

#include <atomic>

namespace ikafssn {

// TempEntry is defined in index/parallel_sort_dispatch.hpp.

// Determine partition for a kmer based on upper bits.
// effective_bits is the total number of significant bits (2*k).
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

    // Pre-compute spaced seed masks (shared across all phases).
    std::vector<uint32_t> seed_masks;
    if (config.t > 0) {
        seed_masks = get_seed_masks(k, config.t,
                         static_cast<TemplateType>(config.template_type));
    }

    const uint32_t tbl_size = table_size(k);
    const int effective_bits = 2 * k;

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
                            config.max_degen_expand);
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

        // Reduce thread-local counts.  Walking each thread-local vector
        // sequentially is O(nthread * tbl_size) on a single thread and
        // becomes a bottleneck at k=13/15 (tbl_size = 4^k = 64M / 1G entries).
        // Snapshot the local pointers, then range-partition the tbl_size
        // dimension across threads — each output cell aggregates across all
        // thread-local snapshots in parallel without touching extra memory
        // (parallel_reduce-style identity-vector splitting would push peak
        // memory to N+log(N) copies of `counts64`, which is intolerable at
        // these table sizes).
        std::vector<const std::vector<uint64_t>*> snaps;
        local_counts.combine_each([&](const std::vector<uint64_t>& lc) {
            snaps.push_back(&lc);
        });

        tbb::parallel_for(
            tbb::blocked_range<uint32_t>(0, tbl_size, 64 * 1024),
            [&](const tbb::blocked_range<uint32_t>& r) {
                for (uint32_t i = r.begin(); i < r.end(); i++) {
                    uint64_t s = 0;
                    for (auto* p : snaps) s += (*p)[i];
                    counts64[i] = s;
                }
            });
    }

    // Convert uint64 -> uint32 with overflow check
    std::vector<uint32_t> counts(tbl_size);
    uint64_t total_occurrences = 0;
    for (uint32_t i = 0; i < tbl_size; i++) {
        if (counts64[i] > UINT32_MAX) {
            logger.error("k-mer %u has count %lu which exceeds uint32_t. "
                         "Use a larger k value.", i, counts64[i]);
            std::remove(ksx_tmp.c_str());
            return false;
        }
        counts[i] = static_cast<uint32_t>(counts64[i]);
        total_occurrences += counts64[i];
    }
    counts64.clear();
    counts64.shrink_to_fit();

    logger.info("Phase 1: total k-mer occurrences = %lu",
                static_cast<unsigned long>(total_occurrences));

    // =========== Determine partition count from memory_limit ===========
    // Per-entry overhead depends on the active sort path: the SIMD path needs
    // a parallel uint64_t key + uint32_t val array on top of the TempEntry
    // buffer (Strategy E in parallel_sort_dispatch.cpp).
    int num_partitions = 1;
    if (total_occurrences > 0) {
        uint64_t entries_limit = config.memory_limit / parallel_sort_entry_overhead();
        while (static_cast<uint64_t>((total_occurrences + num_partitions - 1) / num_partitions)
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

    // posting file starts here (no counts table in format v3)
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

    // Current write positions in the posting file (relative to posting file start)
    uint64_t kix_data_pos = 0;
    uint64_t kpx_data_pos = 0;

    // Per-thread scratch for the parallel per-k-mer encode step below.
    // The buffers grow monotonically across k-mer runs so the dedup +
    // encode pipeline avoids per-call allocation in the hot path.
    struct EncodeScratch {
        std::vector<uint32_t> sid_buf;
        std::vector<uint32_t> abs_pos_buf;
        std::vector<uint32_t> distinct_sid_buf;
        std::vector<uint32_t> occ_count_buf;
        std::vector<uint32_t> seq_delta_buf;
    };

    // Per-k-mer run output: kix / kpx posting list bodies plus the
    // distinct / position counts that the prefix-sum step needs.
    struct RunOut {
        std::vector<uint8_t> kix_blob;
        std::vector<uint8_t> kpx_blob;   // empty when skip_kpx
        uint32_t distinct_count = 0;
        uint32_t position_count = 0;
    };

    // Track v7 totals for the .kix and .kpx headers separately:
    //   total_distinct_postings — sum of distinct seq_ids across all k-mers (.kix)
    //   total_position_count    — sum of intra-sequence position counts (.kpx)
    uint64_t total_distinct_postings = 0;
    uint64_t total_position_count = 0;

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
                            config.max_degen_expand);
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

        // Merge thread-local buffers into a single buffer.  The serial
        // `combine_each(insert)` form was a sustained 1-CPU stretch on
        // large partitions; instead, collect the per-thread vectors,
        // prefix-sum their sizes, allocate the destination once, and
        // copy each thread's slice in parallel.
        std::vector<std::vector<TempEntry>*> local_ptrs;
        local_buffers.combine_each([&](std::vector<TempEntry>& local) {
            local_ptrs.push_back(&local);
        });
        std::vector<std::size_t> chunk_off(local_ptrs.size() + 1, 0);
        for (std::size_t t = 0; t < local_ptrs.size(); ++t) {
            chunk_off[t + 1] = chunk_off[t] + local_ptrs[t]->size();
        }
        std::vector<TempEntry> buffer;
        buffer.resize(chunk_off.back());
        if (!local_ptrs.empty()) {
            tbb::parallel_for(
                tbb::blocked_range<std::size_t>(0, local_ptrs.size(), 1),
                [&](const tbb::blocked_range<std::size_t>& r) {
                    for (std::size_t t = r.begin(); t < r.end(); ++t) {
                        auto& src = *local_ptrs[t];
                        std::move(src.begin(), src.end(),
                                  buffer.begin() + chunk_off[t]);
                        std::vector<TempEntry>().swap(src);
                    }
                });
        }

        if (buffer.empty()) continue;

        // Sort by (kmer_value, seq_id, pos). Dispatched to SIMD path
        // (x86-simd-sort) on AVX2+ x86_64 / TBB parallel_sort otherwise.
        parallel_sort_temp_entries(buffer);

        // Step 1 (sequential): identify k-mer runs in the sorted buffer.
        // Cheap linear scan; the heavy work — dedup + encode_posting_kix /
        // encode_posting_kpx per run — is parallelised below.
        struct KmerRun {
            uint32_t kmer;
            std::size_t begin;
            std::size_t end;
        };
        std::vector<KmerRun> runs;
        {
            std::size_t i = 0;
            while (i < buffer.size()) {
                uint32_t cur_kmer = buffer[i].kmer_value;
                std::size_t j = i;
                while (j < buffer.size() && buffer[j].kmer_value == cur_kmer) {
                    j++;
                }
                // The .kix / .kpx encoder APIs take position_count as u32,
                // so a single (k-mer, partition) slice cannot exceed
                // 2^32 - 1.  In practice this only matters for hypothetical
                // mega-volumes (NCBI BLAST splits at ~4 GB of sequence
                // data, where the dominant k-mer typically has < 100 M
                // occurrences).  Detect the overflow and abort with a
                // clear error rather than silently truncating and
                // corrupting downstream buffers.
                if ((j - i) > UINT32_MAX) {
                    logger.error("k-mer %u has %zu positions in this partition, "
                                 "exceeding uint32_t.  Reduce -memory_limit to "
                                 "force more partitions, or split the BLAST DB "
                                 "into smaller volumes.",
                                 cur_kmer, j - i);
                    std::fclose(kix_fp);
                    if (kpx_fp) std::fclose(kpx_fp);
                    std::remove(ksx_tmp.c_str());
                    std::remove(kix_tmp.c_str());
                    if (!config.skip_kpx) std::remove(kpx_tmp.c_str());
                    return false;
                }
                runs.push_back({cur_kmer, i, j});
                i = j;
            }
        }

        // Step 2 (parallel): per-k-mer dedup + encode into per-run blobs.
        // Each thread reuses an EncodeScratch instance across runs so the
        // dedup and delta-stream allocations don't churn in the hot path.
        // RunOut holds the per-run output buffers; the sequential fwrite
        // loop below walks them in k-mer order.
        std::vector<RunOut> run_outs(runs.size());
        tbb::combinable<EncodeScratch> scratch_ets;

        tbb::parallel_for(
            tbb::blocked_range<std::size_t>(0, runs.size()),
            [&](const tbb::blocked_range<std::size_t>& r) {
                auto& s = scratch_ets.local();
                for (std::size_t ri = r.begin(); ri < r.end(); ri++) {
                    const KmerRun& run = runs[ri];
                    const uint32_t position_count =
                        static_cast<uint32_t>(run.end - run.begin);

                    // Materialise per-k-mer sid / pos arrays (sorted by
                    // (seq_id, pos)).  Both feed the dedup + encode pipeline.
                    s.sid_buf.resize(position_count);
                    s.abs_pos_buf.resize(position_count);
                    for (std::size_t e = 0; e < position_count; e++) {
                        s.sid_buf[e]     = buffer[run.begin + e].seq_id;
                        s.abs_pos_buf[e] = buffer[run.begin + e].pos;
                    }

                    // SIMD dedup → distinct_sid + per-seq_id occurrence count.
                    s.distinct_sid_buf.resize(position_count);
                    s.occ_count_buf.resize(position_count);
                    const uint32_t distinct_count = seq_id_dedup::dedup_seq_ids(
                        s.sid_buf.data(), position_count,
                        s.distinct_sid_buf.data(), s.occ_count_buf.data());
                    s.distinct_sid_buf.resize(distinct_count);
                    s.occ_count_buf.resize(distinct_count);

                    // .kix: distinct seq_id delta-stream → FastPFor.
                    s.seq_delta_buf.resize(distinct_count);
                    if (distinct_count > 0) {
                        s.seq_delta_buf[0] = s.distinct_sid_buf[0];
                        for (uint32_t e = 1; e < distinct_count; e++) {
                            s.seq_delta_buf[e] =
                                s.distinct_sid_buf[e] - s.distinct_sid_buf[e - 1];
                        }
                    }
                    pfd::encode_posting_kix(s.seq_delta_buf.data(),
                                            distinct_count,
                                            run_outs[ri].kix_blob);

                    // .kpx (skip if mode 1).
                    if (!config.skip_kpx) {
                        pfd::encode_posting_kpx(s.distinct_sid_buf.data(),
                                                s.occ_count_buf.data(),
                                                distinct_count,
                                                s.abs_pos_buf.data(),
                                                position_count,
                                                config.freq_threshold_part,
                                                run_outs[ri].kpx_blob);
                    }

                    run_outs[ri].distinct_count = distinct_count;
                    run_outs[ri].position_count = position_count;
                }
            });

        // Step 3 (sequential): record dictionary offsets, accumulate header
        // totals, and stream the per-run blobs to disk in k-mer order.
        // Runs are already sorted by k-mer because `runs` was built from a
        // (kmer, seq_id, pos)-sorted buffer.
        for (std::size_t ri = 0; ri < runs.size(); ri++) {
            const uint32_t cur_kmer = runs[ri].kmer;
            auto& out = run_outs[ri];

            kix_offsets[cur_kmer] = kix_data_pos;
            if (!config.skip_kpx) kpx_offsets[cur_kmer] = kpx_data_pos;
            total_distinct_postings += out.distinct_count;
            total_position_count    += out.position_count;

            std::fwrite(out.kix_blob.data(), 1, out.kix_blob.size(), kix_fp);
            kix_data_pos += out.kix_blob.size();
            std::vector<uint8_t>().swap(out.kix_blob);

            if (!config.skip_kpx) {
                std::fwrite(out.kpx_blob.data(), 1, out.kpx_blob.size(), kpx_fp);
                kpx_data_pos += out.kpx_blob.size();
                std::vector<uint8_t>().swap(out.kpx_blob);
            }
        }

        logger.debug("  Partition %d: %lu entries written", p + 1,
                     static_cast<unsigned long>(buffer.size()));
    }

    // Forward-fill kix_offsets: empty k-mers get the same offset as the next
    // non-empty k-mer (or the sentinel). This ensures offsets[i+1]-offsets[i]==0
    // for empty k-mers and keeps the array monotonically non-decreasing so
    // it can be Elias-Fano encoded (Phase 7a).
    {
        uint64_t fill = kix_data_pos; // sentinel value for trailing empties
        for (int32_t i = static_cast<int32_t>(tbl_size) - 1; i >= 0; i--) {
            if (counts[i] > 0) {
                fill = kix_offsets[i];
            } else {
                kix_offsets[i] = fill;
            }
        }
    }

    // Phase 7e: same treatment for kpx_offsets — needed so the EF
    // encoder sees a monotonically non-decreasing input (the v8 raw
    // u32/u64 dictionary tolerated the leading zeros because pos_offset
    // was only ever queried after kix_len > 0 gating, but EF encoding
    // requires monotonic input).
    if (!config.skip_kpx) {
        uint64_t fill = kpx_data_pos;
        for (int32_t i = static_cast<int32_t>(tbl_size) - 1; i >= 0; i--) {
            if (counts[i] > 0) {
                fill = kpx_offsets[i];
            } else {
                kpx_offsets[i] = fill;
            }
        }
    }

    // =========== Phase 4: Finalize ===========
    logger.info("Phase 4: finalizing (distinct_postings=%lu, position_count=%lu)...",
                static_cast<unsigned long>(total_distinct_postings),
                static_cast<unsigned long>(total_position_count));

    // Set sentinel offset
    kix_offsets[tbl_size] = kix_data_pos;

    // Phase 7a/7e: both .kix and .kpx use Elias-Fano; the offset-width
    // selection branch is gone.  Files are still rewritten at finalize
    // because the in-progress placeholders use raw u64 offsets and need
    // to be replaced with the EF blob.

    // Close the in-progress kix/kpx files; we'll rewrite them below
    // with the EF dictionary.  Both must be flushed before we reopen
    // them for read, otherwise stdio's userspace write buffer can
    // shadow the just-written tail bytes of the last few k-mers.
    std::fclose(kix_fp);
    if (kpx_fp) std::fclose(kpx_fp);

    // Rewrite .kix with correct offset width
    {
        // Read the posting file from the temp file (skip header + uint64 offsets)
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
        kix_hdr.total_distinct_postings = total_distinct_postings;
        // Phase 7a: KIX_FLAG_OFFSET32 reserved (Elias-Fano dictionary follows).
        kix_hdr.flags = KIX_FLAG_HAS_KSX;
        kix_hdr.volume_index = volume_index;
        kix_hdr.total_volumes = total_volumes;
        size_t name_len = std::min(db_name.size(), size_t(32));
        kix_hdr.db_len = static_cast<uint16_t>(name_len);
        std::memcpy(kix_hdr.db, db_name.c_str(), name_len);
        kix_hdr.t = config.t;
        kix_hdr.template_type = config.template_type;
        // Reserved codec-extension area (zero in v5; codec follows
        // format_version since Phase 5g-1).
        kix_hdr.codec_id              = 0;
        kix_hdr.codec_version         = 0;
        kix_hdr.block_size            = 0;
        kix_hdr.tail_codec            = 0;
        kix_hdr.exception_codec_flags = 0;
        std::fwrite(&kix_hdr, sizeof(kix_hdr), 1, wr);

        if (!write_kix_dictionary_ef(wr, kix_offsets.data(), tbl_size,
                                     kix_data_pos)) {
            logger.error("Failed to write EF dictionary to %s", kix_tmp.c_str());
            std::fclose(wr);
            return false;
        }

        if (!posting_blob.empty()) {
            std::fwrite(posting_blob.data(), 1, posting_blob.size(), wr);
        }
        std::fclose(wr);
    }

    // Rewrite .kpx with correct offset width (skip if mode 1)
    if (!config.skip_kpx) {
        // kpx posting file start: header + uint64 offsets
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
        kpx_hdr.total_position_count = total_position_count;
        // Phase 7e: pos_offsets is Elias-Fano; offset_type takes the EF
        // sentinel byte (0xFF) and is no longer consulted at read time.
        kpx_hdr.offset_type = 0xFF;
        // Reserved codec-extension area (zero in v5).
        kpx_hdr.codec_id      = 0;
        kpx_hdr.codec_version = 0;
        kpx_hdr.block_size    = 0;
        kpx_hdr.tail_codec    = 0;
        std::fwrite(&kpx_hdr, sizeof(kpx_hdr), 1, wr);

        if (!write_kpx_dictionary_ef(wr, kpx_offsets.data(), tbl_size,
                                     kpx_data_pos)) {
            logger.error("Failed to write EF dictionary to %s", kpx_tmp.c_str());
            std::fclose(wr);
            return false;
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
