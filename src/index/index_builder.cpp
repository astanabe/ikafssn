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
#include "index/fragment_splitter.hpp"
#include "util/logger.hpp"
#include "util/progress.hpp"

#include <cassert>
#include <cerrno>
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

// A short `fwrite` and a failing `fclose` (the final flush) are how a full
// disk shows up; neither is visible unless the return value is inspected.
static bool write_checked(const void* data, size_t size, size_t count,
                          FILE* fp, const std::string& path,
                          const Logger& logger) {
    if (size == 0 || count == 0) return true;
    if (std::fwrite(data, size, count, fp) != count) {
        logger.error("Failed to write %lu byte(s) to %s: %s",
                     static_cast<unsigned long>(size * count),
                     path.c_str(), std::strerror(errno));
        return false;
    }
    return true;
}

static bool close_checked(FILE* fp, const std::string& path,
                          const Logger& logger) {
    if (std::fclose(fp) != 0) {
        logger.error("Failed to close %s: %s", path.c_str(), std::strerror(errno));
        return false;
    }
    return true;
}

bool build_metadata(BlastDbReader& db,
                    const IndexBuilderConfig& config,
                    const std::string& output_prefix,
                    VolumeMetadata& out,
                    const Logger& logger) {

    const uint32_t blast_num_seqs = db.num_sequences();
    const bool splitting_enabled = (config.min_length_split > 0);

    logger.info("Building metadata: blast-db sequences=%u, "
                "min_seq_length=%u, min_length_split=%u, overlap_length=%u",
                blast_num_seqs,
                config.min_seq_length,
                config.min_length_split,
                config.overlap_length);

    std::string ksx_tmp = output_prefix + ".ksx.tmp";

    // BLAST DB OIDs shorter than min_seq_length are skipped; the rest
    // become parents.  Each parent is split into fragments by
    // fragment_splitter::split() (cutting at ncbi4na==0xF runs so
    // fragments never straddle long ambiguous regions); with
    // min_length_split == 0 every parent emits one whole-parent
    // fragment.  The per-fragment mapping returned via VolumeMetadata
    // lets the postings pass slice each parent's k-mer stream into
    // fragment-relative coordinates.
    logger.info("Collecting metadata%s...",
                splitting_enabled ? " and splitting parents into fragments" : "");
    out = VolumeMetadata{};
    out.seq_id_to_blast_oid.reserve(blast_num_seqs);
    out.seq_id_to_frag_start.reserve(blast_num_seqs);
    out.seq_id_to_frag_end.reserve(blast_num_seqs);

    // The parallel pass below uses only lock-free CSeqDB calls
    // (`seq_length`, `get_raw_sequence`) so it scales freely on a
    // single volume.  `get_accession` (= `CSeqDBImpl::GetSeqIDs`) is
    // deliberately excluded: it acquires the global
    // `CSeqDBAtlas::m_Lock` and walks a non-thread-safe defline
    // cache, so calling it here would serialise the whole loop.
    // It runs instead in the serial merge below, where a single
    // thread holds the lock and touches the cache.
    struct PerOid {
        uint32_t slen = 0;
        bool kept = false;
        std::vector<fragment_splitter::Fragment> frags;
    };
    std::vector<PerOid> per_oid(blast_num_seqs);

    Progress prog("Metadata", blast_num_seqs, config.verbose);
    std::atomic<uint32_t> progress_done{0};
    tbb::parallel_for(
        tbb::blocked_range<uint32_t>(0, blast_num_seqs, 256),
        [&](const tbb::blocked_range<uint32_t>& range) {
            for (uint32_t oid = range.begin(); oid < range.end(); oid++) {
                PerOid& po = per_oid[oid];
                po.slen = db.seq_length(oid);
                if (po.slen < config.min_seq_length) {
                    po.kept = false;
                    continue;
                }
                po.kept = true;
                if (!splitting_enabled) {
                    po.frags.push_back({1u, po.slen});
                } else {
                    auto raw = db.get_raw_sequence(oid);
                    auto ambig = AmbiguityParser::parse(raw.ambig_data,
                                                       raw.ambig_bytes);
                    db.ret_raw_sequence(raw);
                    po.frags = fragment_splitter::split(po.slen, ambig,
                                                        config.min_length_split,
                                                        config.overlap_length);
                }
            }
            uint32_t done = progress_done.fetch_add(
                range.end() - range.begin(), std::memory_order_relaxed)
                + (range.end() - range.begin());
            prog.update(done);
        });
    prog.finish();

    KsxWriter ksx;
    uint32_t skipped = 0;
    for (uint32_t blast_oid = 0; blast_oid < blast_num_seqs; blast_oid++) {
        PerOid& po = per_oid[blast_oid];
        if (!po.kept) {
            skipped++;
            continue;
        }
        std::string acc;
        if (!db.get_accession(blast_oid, acc)) {
            logger.error("Cannot read the accession of BLAST DB OID %u.  A "
                         "common cause is running out of file descriptors "
                         "(see `ulimit -n`); the volume may also be "
                         "truncated or unreadable.", blast_oid);
            std::remove(ksx_tmp.c_str());
            return false;
        }
        uint32_t parent_idx = ksx.add_parent(blast_oid, po.slen, acc);
        for (const auto& f : po.frags) {
            if (out.seq_id_to_blast_oid.size() >= UINT32_MAX) {
                logger.error("SeqId space exhausted: more than %u fragments "
                             "would be required for this volume "
                             "(BLAST DB OID %u, parent length %u). "
                             "Increase -min_length_split or split the BLAST "
                             "DB into smaller volumes.",
                             UINT32_MAX, blast_oid, po.slen);
                std::remove(ksx_tmp.c_str());
                return false;
            }
            uint32_t seq_id = ksx.add_fragment(parent_idx, f.start, f.end);
            (void)seq_id;
            out.seq_id_to_blast_oid.push_back(blast_oid);
            out.seq_id_to_frag_start.push_back(f.start);
            out.seq_id_to_frag_end.push_back(f.end);
        }
        // Drop per-OID frag buffer as we merge to cap peak memory.
        po.frags.clear(); po.frags.shrink_to_fit();
    }

    if (!ksx.write(ksx_tmp)) {
        logger.error("Failed to write %s", ksx_tmp.c_str());
        return false;
    }
    out.num_sequences = ksx.num_sequences();
    out.num_parents = ksx.num_parents();
    logger.info("Wrote %s (%u parents -> %u fragments, "
                "%u skipped < min_seq_length=%u)",
                ksx_tmp.c_str(), out.num_parents, out.num_sequences,
                skipped, config.min_seq_length);
    return true;
}

template <typename KmerInt>
bool build_postings(BlastDbReader& db,
                    const IndexBuilderConfig& config,
                    const std::string& output_prefix,
                    const VolumeMetadata& meta,
                    uint16_t volume_index,
                    uint16_t total_volumes,
                    const std::string& db_name,
                    const Logger& logger) {

    const int k = config.k;

    // File paths (.tmp during construction, renamed to final on success).
    std::string ksx_tmp = output_prefix + ".ksx.tmp";
    std::string kix_tmp = output_prefix + ".kix.tmp";
    std::string kpx_tmp;
    if (!config.skip_kpx) {
        kpx_tmp = output_prefix + ".kpx.tmp";
    }

    const std::vector<uint32_t>& seq_id_to_blast_oid   = meta.seq_id_to_blast_oid;
    const std::vector<uint32_t>& seq_id_to_frag_start  = meta.seq_id_to_frag_start;
    const std::vector<uint32_t>& seq_id_to_frag_end    = meta.seq_id_to_frag_end;

    const uint32_t num_seqs = static_cast<uint32_t>(seq_id_to_blast_oid.size());

    // Pre-compute spaced seed masks (shared across all passes).
    std::vector<uint32_t> seed_masks;
    if (config.t > 0) {
        seed_masks = get_seed_masks(k, config.t,
                         static_cast<TemplateType>(config.template_type));
    }

    const uint32_t tbl_size = table_size(k);
    const int effective_bits = 2 * k;

    // Group fragments by parent so each parent's k-mer stream is
    // scanned once and split into fragment coordinates inline.  All
    // fragments of one parent occupy a contiguous seq_id range with
    // frag_start in ascending order (see build_metadata).
    struct ParentGroup {
        uint32_t blast_oid;
        uint32_t seq_id_lo;  // inclusive
        uint32_t seq_id_hi;  // exclusive
    };
    std::vector<ParentGroup> groups;
    groups.reserve(meta.num_parents);
    for (uint32_t s = 0; s < num_seqs; ) {
        const uint32_t oid = seq_id_to_blast_oid[s];
        uint32_t e = s + 1;
        while (e < num_seqs && seq_id_to_blast_oid[e] == oid) {
            assert(seq_id_to_frag_start[e] >= seq_id_to_frag_start[e - 1]);
            ++e;
        }
        groups.push_back({oid, s, e});
        s = e;
    }
    const uint32_t num_groups = static_cast<uint32_t>(groups.size());

    // Per-fragment parent-relative k-mer pos bounds (0-based leftmost
    // base of the k-mer window).  Fragments shorter than k host no
    // k-mer; for those we set frag_pos_lo to UINT32_MAX so the range
    // check fails unconditionally.
    std::vector<uint32_t> frag_pos_lo(num_seqs);
    std::vector<uint32_t> frag_pos_hi(num_seqs);
    for (uint32_t i = 0; i < num_seqs; i++) {
        const uint32_t fs = seq_id_to_frag_start[i];
        const uint32_t fe = seq_id_to_frag_end[i];
        const bool long_enough =
            (fe >= fs + static_cast<uint32_t>(k) - 1u);
        frag_pos_lo[i] = long_enough ? (fs - 1u) : UINT32_MAX;
        frag_pos_hi[i] = long_enough
            ? (fe - static_cast<uint32_t>(k))
            : 0u;
    }

    // Counting pass (TBB parallel).
    logger.info("Counting k-mers (threads=%d)...", config.threads);
    std::vector<uint64_t> counts64(tbl_size, 0);
    {
        // Always use parallel path; TBB respects global_control parallelism
        // (threads==1 degrades gracefully to sequential execution)
        tbb::combinable<std::vector<uint64_t>> local_counts(
            [&tbl_size]() { return std::vector<uint64_t>(tbl_size, 0); });

        tbb::parallel_for(
            tbb::blocked_range<uint32_t>(0, num_groups, 4),
            [&](const tbb::blocked_range<uint32_t>& range) {
                auto& my_counts = local_counts.local();
                PackedKmerScanner<KmerInt> scanner(k);
                for (uint32_t gi = range.begin(); gi < range.end(); gi++) {
                    const ParentGroup& group = groups[gi];
                    auto raw = db.get_raw_sequence(group.blast_oid);
                    auto ambig = AmbiguityParser::parse(
                        raw.ambig_data, raw.ambig_bytes);
                    // Scanner emits pos non-decreasing, so one
                    // forward pointer resolves covering fragment(s) in
                    // O(1) amortised; overlap_length < min_length_split / 2
                    // (CLI validation) caps the cover at two consecutive
                    // fragments.
                    uint32_t frag_idx = group.seq_id_lo;
                    auto emit_for = [&](uint32_t pos, auto&& on_hit) {
                        while (frag_idx < group.seq_id_hi
                               && pos > frag_pos_hi[frag_idx]) ++frag_idx;
                        if (frag_idx >= group.seq_id_hi) return;
                        if (pos >= frag_pos_lo[frag_idx]) on_hit(frag_idx);
                        const uint32_t nxt = frag_idx + 1;
                        if (nxt < group.seq_id_hi
                            && pos >= frag_pos_lo[nxt]
                            && pos <= frag_pos_hi[nxt]) {
                            on_hit(nxt);
                        }
                    };

                    auto normal_cb = [&](uint32_t pos, KmerInt kmer) {
                        emit_for(pos, [&](uint32_t /*sid*/) {
                            my_counts[kmer]++;
                        });
                    };
                    auto ambig_cb = [&](uint32_t pos, KmerInt base_kmer,
                                        const AmbigInfo* infos, int count) {
                        emit_for(pos, [&](uint32_t /*sid*/) {
                            expand_ambig_kmer_multi<KmerInt>(
                                base_kmer, infos, count,
                                [&](KmerInt expanded) {
                                    my_counts[expanded]++;
                                });
                        });
                    };

                    if (config.t > 0) {
                        scanner.scan_spaced(raw.ncbi2na_data,
                            raw.seq_length, ambig,
                            seed_masks, static_cast<int>(config.t),
                            normal_cb, ambig_cb, config.max_degen_expand);
                    } else {
                        scanner.scan(raw.ncbi2na_data, raw.seq_length,
                            ambig, normal_cb, ambig_cb,
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

    logger.info("Total k-mer occurrences = %lu",
                static_cast<unsigned long>(total_occurrences));

    // Determine partition count from memory_limit.  Per-entry overhead
    // depends on the active sort path: the SIMD path needs a parallel
    // uint64_t key + uint32_t val array on top of the TempEntry buffer.
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
        logger.info("Writing postings (partitions=%d, memory_limit=%luG)...",
                    num_partitions,
                    static_cast<unsigned long>(config.memory_limit >> 30));
    else
        logger.info("Writing postings (partitions=%d, memory_limit=%luM)...",
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

    // Close whatever is still open and drop the partial .tmp files.
    auto abort_postings = [&]() {
        if (kix_fp) { std::fclose(kix_fp); kix_fp = nullptr; }
        if (kpx_fp) { std::fclose(kpx_fp); kpx_fp = nullptr; }
        std::remove(ksx_tmp.c_str());
        std::remove(kix_tmp.c_str());
        if (!config.skip_kpx) std::remove(kpx_tmp.c_str());
        return false;
    };

    // Write kix header placeholder (will be overwritten at finalize)
    KixHeader kix_hdr{};
    if (!write_checked(&kix_hdr, sizeof(kix_hdr), 1, kix_fp, kix_tmp, logger))
        return abort_postings();

    // Reserve the dictionary (u64 offsets here; the filter rewrites it as Elias-Fano)
    const uint64_t kix_offsets_pos = sizeof(KixHeader);
    std::vector<uint64_t> kix_offsets(tbl_size + 1, 0);
    if (!write_checked(kix_offsets.data(), sizeof(uint64_t), tbl_size + 1,
                       kix_fp, kix_tmp, logger))
        return abort_postings();

    // posting file starts here
    const uint64_t kix_posting_start = kix_offsets_pos + sizeof(uint64_t) * (tbl_size + 1);

    // Write kpx header placeholder (skip if mode 1)
    KpxHeader kpx_hdr{};
    std::vector<uint64_t> kpx_offsets;
    if (!config.skip_kpx) {
        if (!write_checked(&kpx_hdr, sizeof(kpx_hdr), 1, kpx_fp, kpx_tmp, logger))
            return abort_postings();

        // Reserve the .kpx dictionary (pos_offsets)
        kpx_offsets.resize(tbl_size, 0);
        if (!write_checked(kpx_offsets.data(), sizeof(uint64_t), tbl_size,
                           kpx_fp, kpx_tmp, logger))
            return abort_postings();
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

    // Header totals tracked separately for .kix and .kpx:
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
            tbb::blocked_range<uint32_t>(0, num_groups, 4),
            [&](const tbb::blocked_range<uint32_t>& range) {
                auto& my_buffer = local_buffers.local();
                PackedKmerScanner<KmerInt> scanner(k);
                for (uint32_t gi = range.begin(); gi < range.end(); gi++) {
                    const ParentGroup& group = groups[gi];
                    auto raw = db.get_raw_sequence(group.blast_oid);
                    auto ambig = AmbiguityParser::parse(
                        raw.ambig_data, raw.ambig_bytes);
                    uint32_t frag_idx = group.seq_id_lo;
                    auto emit_for = [&](uint32_t pos, auto&& on_hit) {
                        while (frag_idx < group.seq_id_hi
                               && pos > frag_pos_hi[frag_idx]) ++frag_idx;
                        if (frag_idx >= group.seq_id_hi) return;
                        if (pos >= frag_pos_lo[frag_idx]) on_hit(frag_idx);
                        const uint32_t nxt = frag_idx + 1;
                        if (nxt < group.seq_id_hi
                            && pos >= frag_pos_lo[nxt]
                            && pos <= frag_pos_hi[nxt]) {
                            on_hit(nxt);
                        }
                    };

                    // .kpx positions are fragment-relative; subtract
                    // frag_pos_lo[sid] from the parent-relative pos.
                    auto normal_cb = [&](uint32_t pos, KmerInt kmer) {
                        const uint32_t kval = static_cast<uint32_t>(kmer);
                        if (counts[kval] == 0) return;
                        if (static_cast<int>(partition_of(
                                kval, partition_bits, effective_bits)) != p)
                            return;
                        emit_for(pos, [&](uint32_t sid) {
                            my_buffer.push_back({kval, sid,
                                pos - frag_pos_lo[sid]});
                        });
                    };
                    auto ambig_cb = [&](uint32_t pos, KmerInt base_kmer,
                                        const AmbigInfo* infos, int count) {
                        emit_for(pos, [&](uint32_t sid) {
                            const uint32_t shift = frag_pos_lo[sid];
                            expand_ambig_kmer_multi<KmerInt>(
                                base_kmer, infos, count,
                                [&](KmerInt expanded) {
                                    const uint32_t kval =
                                        static_cast<uint32_t>(expanded);
                                    if (counts[kval] == 0) return;
                                    if (static_cast<int>(partition_of(
                                            kval, partition_bits,
                                            effective_bits)) != p) return;
                                    my_buffer.push_back({kval, sid,
                                        pos - shift});
                                });
                        });
                    };

                    if (config.t > 0) {
                        scanner.scan_spaced(raw.ncbi2na_data,
                            raw.seq_length, ambig,
                            seed_masks, static_cast<int>(config.t),
                            normal_cb, ambig_cb, config.max_degen_expand);
                    } else {
                        scanner.scan(raw.ncbi2na_data, raw.seq_length,
                            ambig, normal_cb, ambig_cb,
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
                                 100.0 * done / num_groups, done, num_groups,
                                 static_cast<long>(elapsed));
                    std::fflush(stderr);
                }
            });

        if (config.verbose) {
            auto now = std::chrono::steady_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(
                now - progress_start).count();
            std::fprintf(stderr, "\r  Partition scan: done (%u parents, %lds)\n",
                         num_groups, static_cast<long>(elapsed));
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
                // The .kix / .kpx encoder APIs take position_count as u32, so
                // abort if a single (k-mer, partition) slice exceeds 2^32 - 1
                // rather than silently truncating and corrupting downstream
                // buffers.
                if ((j - i) > UINT32_MAX) {
                    logger.error("k-mer %u has %zu positions in this partition, "
                                 "exceeding uint32_t.  Reduce -memory_limit to "
                                 "force more partitions, or split the BLAST DB "
                                 "into smaller volumes.",
                                 cur_kmer, j - i);
                    return abort_postings();
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

            if (!write_checked(out.kix_blob.data(), 1, out.kix_blob.size(),
                               kix_fp, kix_tmp, logger))
                return abort_postings();
            kix_data_pos += out.kix_blob.size();
            std::vector<uint8_t>().swap(out.kix_blob);

            if (!config.skip_kpx) {
                if (!write_checked(out.kpx_blob.data(), 1, out.kpx_blob.size(),
                                   kpx_fp, kpx_tmp, logger))
                    return abort_postings();
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
    // it can be Elias-Fano encoded.
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

    // Same forward-fill for kpx_offsets — the EF encoder requires a
    // monotonically non-decreasing input.
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

    // Finalize: rewrite .kix / .kpx with the Elias-Fano dictionary in
    // place of the raw u64 offset placeholders.
    logger.info("Finalizing (distinct_postings=%lu, position_count=%lu)...",
                static_cast<unsigned long>(total_distinct_postings),
                static_cast<unsigned long>(total_position_count));

    // Set sentinel offset
    kix_offsets[tbl_size] = kix_data_pos;

    // Close the in-progress kix/kpx files; we'll rewrite them below
    // with the EF dictionary.  Both must be flushed before we reopen
    // them for read, otherwise stdio's userspace write buffer can
    // shadow the just-written tail bytes of the last few k-mers.
    {
        FILE* fp = kix_fp;
        kix_fp = nullptr;
        if (!close_checked(fp, kix_tmp, logger)) return abort_postings();
    }
    if (kpx_fp) {
        FILE* fp = kpx_fp;
        kpx_fp = nullptr;
        if (!close_checked(fp, kpx_tmp, logger)) return abort_postings();
    }

    // Rewrite .kix with correct offset width
    {
        // Read the posting file from the temp file (skip header + uint64 offsets)
        FILE* rd = std::fopen(kix_tmp.c_str(), "rb");
        if (!rd) {
            logger.error("Cannot open %s for reading: %s",
                         kix_tmp.c_str(), std::strerror(errno));
            return abort_postings();
        }
        if (std::fseek(rd, static_cast<long>(kix_posting_start), SEEK_SET) != 0) {
            logger.error("Cannot seek to the posting file in %s: %s",
                         kix_tmp.c_str(), std::strerror(errno));
            std::fclose(rd);
            return abort_postings();
        }
        std::vector<uint8_t> posting_file(kix_data_pos);
        if (kix_data_pos > 0 &&
            std::fread(posting_file.data(), 1, kix_data_pos, rd) != kix_data_pos) {
            logger.error("Failed to read the posting file back from %s", kix_tmp.c_str());
            std::fclose(rd);
            return abort_postings();
        }
        std::fclose(rd);

        // Write the final kix file
        FILE* wr = std::fopen(kix_tmp.c_str(), "wb");
        if (!wr) {
            logger.error("Cannot open %s for writing: %s",
                         kix_tmp.c_str(), std::strerror(errno));
            return abort_postings();
        }

        std::memcpy(kix_hdr.magic, KIX_MAGIC, sizeof(KIX_MAGIC));
        kix_hdr.format_version = KIX_FORMAT_VERSION;
        kix_hdr.k = static_cast<uint8_t>(k);
        kix_hdr.kmer_type = kmer_type_for(k, config.t);
        kix_hdr.num_sequences = num_seqs;
        kix_hdr.total_distinct_postings = total_distinct_postings;
        // KIX_FLAG_OFFSET32 reserved; Elias-Fano dictionary follows.
        kix_hdr.flags = KIX_FLAG_HAS_KSX;
        kix_hdr.volume_index = volume_index;
        kix_hdr.total_volumes = total_volumes;
        size_t name_len = std::min(db_name.size(), size_t(32));
        kix_hdr.db_len = static_cast<uint16_t>(name_len);
        std::memcpy(kix_hdr.db, db_name.c_str(), name_len);
        kix_hdr.t = config.t;
        kix_hdr.template_type = config.template_type;
        // Reserved codec-extension area (codec follows format_version).
        kix_hdr.codec_id              = 0;
        kix_hdr.codec_version         = 0;
        kix_hdr.block_size            = 0;
        kix_hdr.tail_codec            = 0;
        kix_hdr.exception_codec_flags = 0;
        if (!write_checked(&kix_hdr, sizeof(kix_hdr), 1, wr, kix_tmp, logger)) {
            std::fclose(wr);
            return abort_postings();
        }

        if (!write_kix_dictionary_ef(wr, kix_offsets.data(), tbl_size,
                                     kix_data_pos)) {
            logger.error("Failed to write EF dictionary to %s", kix_tmp.c_str());
            std::fclose(wr);
            return abort_postings();
        }

        if (!write_checked(posting_file.data(), 1, posting_file.size(),
                           wr, kix_tmp, logger)) {
            std::fclose(wr);
            return abort_postings();
        }
        if (!close_checked(wr, kix_tmp, logger)) return abort_postings();
    }

    // Rewrite .kpx with correct offset width (skip if mode 1)
    if (!config.skip_kpx) {
        // kpx posting file start: header + uint64 offsets
        const uint64_t kpx_posting_start_pos = sizeof(KpxHeader) + sizeof(uint64_t) * tbl_size;

        FILE* rd = std::fopen(kpx_tmp.c_str(), "rb");
        if (!rd) {
            logger.error("Cannot open %s for reading: %s",
                         kpx_tmp.c_str(), std::strerror(errno));
            return abort_postings();
        }
        if (std::fseek(rd, static_cast<long>(kpx_posting_start_pos), SEEK_SET) != 0) {
            logger.error("Cannot seek to the posting file in %s: %s",
                         kpx_tmp.c_str(), std::strerror(errno));
            std::fclose(rd);
            return abort_postings();
        }
        std::vector<uint8_t> posting_file(kpx_data_pos);
        if (kpx_data_pos > 0 &&
            std::fread(posting_file.data(), 1, kpx_data_pos, rd) != kpx_data_pos) {
            logger.error("Failed to read the posting file back from %s", kpx_tmp.c_str());
            std::fclose(rd);
            return abort_postings();
        }
        std::fclose(rd);

        FILE* wr = std::fopen(kpx_tmp.c_str(), "wb");
        if (!wr) {
            logger.error("Cannot open %s for writing: %s",
                         kpx_tmp.c_str(), std::strerror(errno));
            return abort_postings();
        }

        std::memcpy(kpx_hdr.magic, KPX_MAGIC, sizeof(KPX_MAGIC));
        kpx_hdr.format_version = KPX_FORMAT_VERSION;
        kpx_hdr.k = static_cast<uint8_t>(k);
        kpx_hdr.t = config.t;
        kpx_hdr.template_type = config.template_type;
        kpx_hdr.total_position_count = total_position_count;
        // pos_offsets is Elias-Fano; offset_type carries the EF
        // sentinel byte (0xFF) and is not consulted at read time.
        kpx_hdr.offset_type = 0xFF;
        // Reserved codec-extension area.
        kpx_hdr.codec_id      = 0;
        kpx_hdr.codec_version = 0;
        kpx_hdr.block_size    = 0;
        kpx_hdr.tail_codec    = 0;
        if (!write_checked(&kpx_hdr, sizeof(kpx_hdr), 1, wr, kpx_tmp, logger)) {
            std::fclose(wr);
            return abort_postings();
        }

        if (!write_kpx_dictionary_ef(wr, kpx_offsets.data(), tbl_size,
                                     kpx_data_pos)) {
            logger.error("Failed to write EF dictionary to %s", kpx_tmp.c_str());
            std::fclose(wr);
            return abort_postings();
        }

        if (!write_checked(posting_file.data(), 1, posting_file.size(),
                           wr, kpx_tmp, logger)) {
            std::fclose(wr);
            return abort_postings();
        }
        if (!close_checked(wr, kpx_tmp, logger)) return abort_postings();
    }

    logger.info("Postings built: %s (.kix.tmp%s, .ksx.tmp)", output_prefix.c_str(),
                config.skip_kpx ? "" : ", .kpx.tmp");
    return true;
}

// Explicit template instantiations
template bool build_postings<uint16_t>(BlastDbReader&, const IndexBuilderConfig&,
    const std::string&, const VolumeMetadata&, uint16_t, uint16_t,
    const std::string&, const Logger&);
template bool build_postings<uint32_t>(BlastDbReader&, const IndexBuilderConfig&,
    const std::string&, const VolumeMetadata&, uint16_t, uint16_t,
    const std::string&, const Logger&);

template <typename KmerInt>
bool build_index(BlastDbReader& db,
                 const IndexBuilderConfig& config,
                 const std::string& output_prefix,
                 uint16_t volume_index,
                 uint16_t total_volumes,
                 const std::string& db_name,
                 const Logger& logger) {
    VolumeMetadata meta;
    if (!build_metadata(db, config, output_prefix, meta, logger)) {
        return false;
    }
    if (!build_postings<KmerInt>(db, config, output_prefix, meta,
                                 volume_index, total_volumes, db_name, logger)) {
        return false;
    }
    if (config.keep_tmp) return true;

    auto rename_one = [&](const char* ext) {
        std::string from = output_prefix + ext + ".tmp";
        std::string to   = output_prefix + ext;
        if (std::rename(from.c_str(), to.c_str()) != 0) {
            logger.error("Failed to rename %s -> %s", from.c_str(), to.c_str());
            return false;
        }
        return true;
    };
    if (!rename_one(".ksx")) return false;
    if (!rename_one(".kix")) return false;
    if (!config.skip_kpx) {
        if (!rename_one(".kpx")) return false;
    }
    return true;
}

template bool build_index<uint16_t>(BlastDbReader&, const IndexBuilderConfig&,
    const std::string&, uint16_t, uint16_t, const std::string&, const Logger&);
template bool build_index<uint32_t>(BlastDbReader&, const IndexBuilderConfig&,
    const std::string&, uint16_t, uint16_t, const std::string&, const Logger&);

} // namespace ikafssn
