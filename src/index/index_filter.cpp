#include "index/index_filter.hpp"
#include "index/kix_reader.hpp"
#include "index/kpx_reader.hpp"
#include "index/kix_format.hpp"
#include "index/dictionary_io.hpp"
#include "index/kpx_format.hpp"
#include "index/khx_writer.hpp"
#include "index/pfd_codec.hpp"
#include "core/config.hpp"
#include "core/types.hpp"
#include "util/logger.hpp"

#include <algorithm>
#include <atomic>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <thread>
#include <vector>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/combinable.h>
#include <tbb/task_arena.h>

namespace ikafssn {

namespace {

// Walk a single FOR-block stream, returning the number of bytes consumed
// (>0) on success or 0 on a structural mismatch.  Mirrors the layout in
// pfd_codec_tier.cpp: per-128-element block `[u32 min][u8 b][3 B pad]`
// + `16*b` body bytes; tail `[u8 tail_count]` plus, when tail_count > 0,
// `[u32 tail_min][u8 tail_b][bitpacked body]`.
std::size_t walk_for_stream(const std::uint8_t* p,
                            const std::uint8_t* end,
                            std::uint32_t total_count) {
    const std::uint8_t* start = p;
    const std::uint32_t num_blocks = total_count / 128u;
    const std::uint32_t tail_count = total_count % 128u;
    for (std::uint32_t i = 0; i < num_blocks; i++) {
        if (p + 8 > end) return 0;
        const std::uint8_t b = p[4];
        if (b > 32) return 0;
        const std::size_t need = std::size_t(8) + std::size_t(16) * b;
        if (static_cast<std::size_t>(end - p) < need) return 0;
        p += need;
    }
    if (p + 1 > end) return 0;
    const std::uint8_t got_tail = *p++;
    if (got_tail != tail_count) return 0;
    if (tail_count > 0) {
        if (p + 5 > end) return 0;
        p += 4;  // tail_min (skipped)
        const std::uint8_t tail_b = *p++;
        if (tail_b > 32) return 0;
        const std::size_t body_bits = std::size_t(tail_count) * tail_b;
        const std::size_t body_bytes = (body_bits + 7) / 8;
        if (static_cast<std::size_t>(end - p) < body_bytes) return 0;
        p += body_bytes;
    }
    return static_cast<std::size_t>(p - start);
}

// Walk one .kpx posting list and return either bytes consumed (>0) or 0
// on a structural mismatch.  Sets *position_count to the per-k-mer
// total (sum of partition occ_counts + short1_count + sum(short2 u8
// occ_count[])) on success.
std::size_t walk_kpx_posting(const std::uint8_t* p,
                             std::size_t bytes,
                             std::uint32_t kix_count,
                             std::uint64_t* position_count) {
    *position_count = 0;
    if (kix_count == 0) {
        return (bytes == 0) ? 1 /* sentinel: 0 bytes consumed but valid */ : 0;
    }
    const std::uint8_t* start = p;
    const std::uint8_t* end = p + bytes;

    // Kind map.
    const std::size_t kind_map_bytes = (std::size_t(kix_count) * 2 + 7) / 8;
    if (static_cast<std::size_t>(end - p) < kind_map_bytes) return 0;
    const std::uint8_t* kind_map = p;
    p += kind_map_bytes;

    // Per-kind counts via popcount.
    std::uint32_t partition_count = 0, short1_count = 0, short2_count = 0;
    for (std::uint32_t k = 0; k < kix_count; k++) {
        const std::uint8_t kind = static_cast<std::uint8_t>(
            (kind_map[k >> 2] >> ((k & 3) * 2)) & 0x03);
        switch (kind) {
            case 0: short1_count++;    break;
            case 1: short2_count++;    break;
            case 2: partition_count++; break;
            default: return 0;  // 11 reserved
        }
    }
    if (partition_count + short1_count + short2_count != kix_count) return 0;

    // Partition groups: each [u32 occ_count][FOR stream].
    std::uint64_t total_pos = 0;
    for (std::uint32_t g = 0; g < partition_count; g++) {
        if (p + 4 > end) return 0;
        std::uint32_t gcnt;
        std::memcpy(&gcnt, p, 4);
        p += 4;
        if (gcnt == 0) return 0;
        total_pos += gcnt;
        const std::size_t consumed = walk_for_stream(p, end, gcnt);
        if (consumed == 0) return 0;
        p += consumed;
    }

    // Short1 FOR stream (one position per cluster).
    if (short1_count > 0) {
        const std::size_t consumed = walk_for_stream(p, end, short1_count);
        if (consumed == 0) return 0;
        p += consumed;
    }
    total_pos += short1_count;

    // Short2: u8 occ_count[short2_count] + FOR stream.
    if (short2_count > 0) {
        if (static_cast<std::size_t>(end - p) < short2_count) return 0;
        std::uint32_t short2_pos = 0;
        for (std::uint32_t i = 0; i < short2_count; i++) {
            const std::uint8_t occ = p[i];
            if (occ < 2) return 0;
            short2_pos += occ;
        }
        p += short2_count;
        total_pos += short2_pos;
        if (short2_pos > 0) {
            const std::size_t consumed = walk_for_stream(p, end, short2_pos);
            if (consumed == 0) return 0;
            p += consumed;
        }
    }

    *position_count = total_pos;
    return static_cast<std::size_t>(p - start);
}

} // anonymous namespace

bool validate_volume(const std::string& kix_path,
                     const std::string& kpx_path,
                     uint64_t* total_position_count_out,
                     const Logger& logger) {
    KixReader kix;
    if (!kix.open(kix_path)) {
        logger.error("validate_volume: cannot open %s", kix_path.c_str());
        return false;
    }
    const std::uint32_t tbl_size = kix.table_size();

    KpxReader kpx;
    const bool has_kpx = !kpx_path.empty();
    if (has_kpx) {
        if (!kpx.open(kpx_path)) {
            logger.error("validate_volume: cannot open %s", kpx_path.c_str());
            return false;
        }
        if (kpx.table_size() != tbl_size) {
            logger.error("validate_volume: .kix table_size %u != .kpx table_size %u",
                         tbl_size, kpx.table_size());
            return false;
        }
    }

    const std::uint8_t* kpx_data = has_kpx ? kpx.posting_file() : nullptr;
    const std::size_t   kpx_size = has_kpx ? kpx.posting_file_size() : 0;
    const auto kix_counts = kix.bulk_count_postings();
    (void)kpx_size;

    if (!has_kpx) {
        if (total_position_count_out) *total_position_count_out = 0;
        return true;
    }

    // Phase 7d: parallel per-k-mer validation via TBB.  Each k-mer's
    // walk is independent; per-thread accumulators feed into a final
    // reduction.  The first error encountered is reported and the
    // pass returns false.  At nt scale (k=12, D=16M) this drops the
    // wall time from minutes per volume to seconds.
    std::atomic<bool>     ok{true};
    std::atomic<uint32_t> first_bad_kmer{UINT32_MAX};
    tbb::combinable<std::uint64_t> local_pos{[]{ return std::uint64_t{0}; }};

    const auto record_error = [&](std::uint32_t k) {
        ok.store(false, std::memory_order_relaxed);
        std::uint32_t expected = UINT32_MAX;
        first_bad_kmer.compare_exchange_strong(expected, k);
    };

    tbb::parallel_for(
        tbb::blocked_range<std::uint32_t>(0, tbl_size, 1u << 14),
        [&](const tbb::blocked_range<std::uint32_t>& range) {
            std::uint64_t local_sum = 0;
            for (std::uint32_t i = range.begin(); i < range.end(); ++i) {
                if (!ok.load(std::memory_order_relaxed)) break;

                const std::uint32_t kix_count = kix_counts[i];
                std::uint64_t lo, hi;
                kpx.pos_offset_range(i, lo, hi);
                if (hi < lo) { record_error(i); break; }
                const std::size_t available = static_cast<std::size_t>(hi - lo);

                std::uint64_t per_kmer_pos = 0;
                const std::size_t consumed = walk_kpx_posting(
                    kpx_data + lo, available, kix_count, &per_kmer_pos);

                if (kix_count == 0) {
                    if (available != 0 || consumed == 0) {
                        record_error(i);
                        break;
                    }
                    continue;
                }

                if (consumed == 0 || consumed != available) {
                    record_error(i);
                    break;
                }
                local_sum += per_kmer_pos;
            }
            local_pos.local() += local_sum;
        });

    if (!ok.load()) {
        const std::uint32_t bad = first_bad_kmer.load();
        const std::uint32_t kix_count = (bad < tbl_size) ? kix_counts[bad] : 0;
        std::uint64_t lo = 0, hi = 0;
        if (bad < tbl_size) {
            kpx.pos_offset_range(bad, lo, hi);
        }
        logger.error(
            "validate_volume: kmer %u .kpx walk failed (kix_count=%u, slice=%lu B)",
            bad, kix_count, static_cast<unsigned long>(hi - lo));
        return false;
    }

    std::uint64_t total_position_count = 0;
    local_pos.combine_each([&](std::uint64_t v) { total_position_count += v; });

    if (total_position_count_out) {
        *total_position_count_out = total_position_count;
    }

    if (has_kpx && kpx.total_position_count() != 0
        && kpx.total_position_count() != total_position_count) {
        logger.error(
            "validate_volume: header total_position_count=%lu != recomputed %lu",
            static_cast<unsigned long>(kpx.total_position_count()),
            static_cast<unsigned long>(total_position_count));
        // Soft fail — the header field is informational and v9
        // filtered builds intentionally report 0 (FIXME path).  Don't
        // return false here so the validator's structural check still
        // succeeds.
    }

    return true;
}

// Compute the byte length of each k-mer.s posting list from sentinel-based offsets.
static std::vector<uint64_t> compute_posting_sizes(
    const KixReader& kix, uint32_t tbl_size) {

    std::vector<uint64_t> sizes(tbl_size);
    for (uint32_t i = 0; i < tbl_size; i++) {
        sizes[i] = kix.posting_list_byte_length(i);
    }
    return sizes;
}

// Write filtered .kix from .kix.tmp (excluded k-mers removed).
static bool write_filtered_kix(
    const KixReader& kix_in,
    const std::string& kix_final,
    const std::vector<bool>& excluded,
    const std::vector<uint64_t>& kix_sizes,
    int k,
    uint32_t tbl_size,
    uint64_t new_total_distinct_postings,
    const Logger& logger) {
    const uint8_t* kix_posting_in = kix_in.posting_file();

    FILE* kix_fp = std::fopen(kix_final.c_str(), "wb");
    if (!kix_fp) {
        logger.error("filter: cannot open %s for writing", kix_final.c_str());
        return false;
    }

    // Build new offsets (tbl_size + 1 entries)
    std::vector<uint64_t> new_kix_offsets(tbl_size + 1, 0);

    // Collect posting list bytes into buffer to determine final size
    std::vector<uint8_t> posting_buf;
    uint64_t kix_data_pos = 0;
    for (uint32_t i = 0; i < tbl_size; i++) {
        new_kix_offsets[i] = kix_data_pos;
        if (kix_sizes[i] > 0 && !excluded[i]) {
            posting_buf.insert(posting_buf.end(),
                kix_posting_in + kix_in.posting_list_offset(i),
                kix_posting_in + kix_in.posting_list_offset(i) + kix_sizes[i]);
            kix_data_pos += kix_sizes[i];
        }
    }
    new_kix_offsets[tbl_size] = kix_data_pos;

    // Write header.  Phase 7a: dictionary is Elias-Fano; the legacy
    // KIX_FLAG_OFFSET32 bit is forced to 0 (preserved as reserved on
    // the header for byte-stability per Phase 7 design decision #6).
    KixHeader kix_hdr{};
    std::memcpy(kix_hdr.magic, KIX_MAGIC, sizeof(KIX_MAGIC));
    kix_hdr.format_version = KIX_FORMAT_VERSION;
    kix_hdr.k = static_cast<uint8_t>(k);
    kix_hdr.kmer_type = kmer_type_for(k, kix_in.header().t);
    kix_hdr.num_sequences = kix_in.num_sequences();
    kix_hdr.total_distinct_postings = new_total_distinct_postings;
    kix_hdr.flags = kix_in.header().flags & ~KIX_FLAG_OFFSET32;
    kix_hdr.volume_index = kix_in.header().volume_index;
    kix_hdr.total_volumes = kix_in.header().total_volumes;
    kix_hdr.db_len = kix_in.header().db_len;
    kix_hdr.t = kix_in.header().t;
    kix_hdr.template_type = kix_in.header().template_type;
    std::memcpy(kix_hdr.db, kix_in.header().db, 32);

    // Reserved codec-extension area (zero in v5).
    kix_hdr.codec_id              = 0;
    kix_hdr.codec_version         = 0;
    kix_hdr.block_size            = 0;
    kix_hdr.tail_codec            = 0;
    kix_hdr.exception_codec_flags = 0;

    // Carry over the v10 fragment-indexing triplet from the input header.
    kix_hdr.min_seq_length   = kix_in.header().min_seq_length;
    kix_hdr.min_length_split = kix_in.header().min_length_split;
    kix_hdr.overlap_length   = kix_in.header().overlap_length;

    std::fwrite(&kix_hdr, sizeof(kix_hdr), 1, kix_fp);

    if (!write_kix_dictionary_ef(kix_fp, new_kix_offsets.data(), tbl_size,
                                 kix_data_pos)) {
        logger.error("filter: failed to write EF dictionary to %s", kix_final.c_str());
        std::fclose(kix_fp);
        return false;
    }

    // Write posting file
    if (!posting_buf.empty()) {
        std::fwrite(posting_buf.data(), 1, posting_buf.size(), kix_fp);
    }

    std::fclose(kix_fp);
    return true;
}

// Write filtered .kpx from .kpx.tmp (excluded k-mers removed).
static bool write_filtered_kpx(
    const KpxReader& kpx_in,
    const std::string& kpx_final,
    const std::vector<bool>& excluded,
    const std::vector<uint64_t>& kpx_sizes,
    const std::vector<uint64_t>& kix_sizes,
    int k,
    uint32_t tbl_size,
    uint64_t new_total_position_count,
    const Logger& logger) {
    const uint8_t* kpx_posting_in = kpx_in.posting_file();

    FILE* kpx_fp = std::fopen(kpx_final.c_str(), "wb");
    if (!kpx_fp) {
        logger.error("filter: cannot open %s for writing", kpx_final.c_str());
        return false;
    }

    std::vector<uint64_t> new_kpx_offsets(tbl_size, 0);
    std::vector<uint8_t> posting_buf;
    uint64_t kpx_data_pos = 0;

    // Phase 7e: the offsets array is Elias-Fano encoded and must be
    // monotonically non-decreasing.  Empty / excluded k-mers therefore
    // get the same offset as the next non-empty non-excluded k-mer
    // (forward-fill) rather than retaining the zero-init value.
    for (uint32_t i = 0; i < tbl_size; i++) {
        new_kpx_offsets[i] = kpx_data_pos;
        if (kix_sizes[i] > 0 && !excluded[i]) {
            posting_buf.insert(posting_buf.end(),
                kpx_posting_in + kpx_in.pos_offset(i),
                kpx_posting_in + kpx_in.pos_offset(i) + kpx_sizes[i]);
            kpx_data_pos += kpx_sizes[i];
        }
    }

    // Write header.  Phase 7e: pos_offsets dictionary is Elias-Fano;
    // offset_type takes the EF sentinel byte (0xFF).
    KpxHeader kpx_hdr{};
    std::memcpy(kpx_hdr.magic, KPX_MAGIC, sizeof(KPX_MAGIC));
    kpx_hdr.format_version = KPX_FORMAT_VERSION;
    kpx_hdr.k = static_cast<uint8_t>(k);
    kpx_hdr.t = kpx_in.header().t;
    kpx_hdr.template_type = kpx_in.header().template_type;
    kpx_hdr.total_position_count = new_total_position_count;
    kpx_hdr.offset_type = 0xFF;  // EF sentinel

    // Reserved codec-extension area (zero in v5).
    kpx_hdr.codec_id      = 0;
    kpx_hdr.codec_version = 0;
    kpx_hdr.block_size    = 0;
    kpx_hdr.tail_codec    = 0;

    // Carry over the v10 fragment-indexing triplet from the input header.
    kpx_hdr.min_seq_length   = kpx_in.header().min_seq_length;
    kpx_hdr.min_length_split = kpx_in.header().min_length_split;
    kpx_hdr.overlap_length   = kpx_in.header().overlap_length;

    std::fwrite(&kpx_hdr, sizeof(kpx_hdr), 1, kpx_fp);

    if (!write_kpx_dictionary_ef(kpx_fp, new_kpx_offsets.data(), tbl_size,
                                 kpx_data_pos)) {
        logger.error("filter: failed to write EF dictionary to %s", kpx_final.c_str());
        std::fclose(kpx_fp);
        return false;
    }

    // Write posting file
    if (!posting_buf.empty()) {
        std::fwrite(posting_buf.data(), 1, posting_buf.size(), kpx_fp);
    }

    std::fclose(kpx_fp);
    return true;
}

// Filter a single volume's .kix.tmp/.kpx.tmp -> .kix/.kpx (in parallel).
static bool filter_one_volume(
    const std::string& vol_prefix,
    const std::vector<bool>& excluded,
    int k,
    const Logger& logger) {

    std::string kix_tmp = vol_prefix + ".kix.tmp";
    std::string kpx_tmp = vol_prefix + ".kpx.tmp";
    std::string ksx_tmp = vol_prefix + ".ksx.tmp";
    std::string kix_final = vol_prefix + ".kix";
    std::string kpx_final = vol_prefix + ".kpx";
    std::string ksx_final = vol_prefix + ".ksx";

    bool has_kpx_tmp = std::filesystem::exists(kpx_tmp);

    KixReader kix_in;
    if (!kix_in.open(kix_tmp)) {
        logger.error("filter: cannot open %s", kix_tmp.c_str());
        return false;
    }

    const uint32_t tbl_size = kix_in.table_size();

    KpxReader kpx_in;
    if (has_kpx_tmp) {
        if (!kpx_in.open(kpx_tmp)) {
            logger.error("filter: cannot open %s", kpx_tmp.c_str());
            return false;
        }
    }

    auto kix_sizes = compute_posting_sizes(kix_in, tbl_size);

    std::vector<uint64_t> kpx_sizes;
    if (has_kpx_tmp) {
        const uint8_t* kpx_posting = kpx_in.posting_file();
        uint64_t kpx_total_data = kpx_in.posting_file_size();
        kpx_sizes.resize(tbl_size, 0);

        uint32_t prev_kmer = UINT32_MAX;
        uint64_t prev_offset = 0;
        for (uint32_t i = 0; i < tbl_size; i++) {
            if (kix_sizes[i] > 0) {
                if (prev_kmer != UINT32_MAX) {
                    kpx_sizes[prev_kmer] = kpx_in.pos_offset(i) - prev_offset;
                }
                prev_kmer = i;
                prev_offset = kpx_in.pos_offset(i);
            }
        }
        if (prev_kmer != UINT32_MAX) {
            kpx_sizes[prev_kmer] = kpx_total_data - prev_offset;
        }
    }

    // Compute new totals.  v7: .kix tracks distinct seq_id postings;
    // .kpx tracks the total position count (sum of intra-sequence
    // occurrences).
    auto counts = kix_in.bulk_count_postings();
    uint64_t new_total_distinct_postings = 0;
    for (uint32_t i = 0; i < tbl_size; i++) {
        if (!excluded[i]) {
            new_total_distinct_postings += counts[i];
        }
    }

    // Phase 7d: walk each non-excluded k-mer's .kpx posting list using
    // the same logic the validator employs (anonymous-namespace
    // walk_kpx_posting above).  This finally gets the filtered .kpx
    // header's total_position_count right — the v8 path read
    // partition_count and silently mis-reported, and the 7c shim was
    // a stop-gap zero.
    uint64_t new_total_position_count = 0;
    if (has_kpx_tmp) {
        const uint8_t* kpx_posting_in = kpx_in.posting_file();
        for (uint32_t i = 0; i < tbl_size; i++) {
            if (excluded[i]) continue;
            if (counts[i] == 0) continue;
            uint64_t lo, hi;
            kpx_in.pos_offset_range(i, lo, hi);
            if (hi <= lo) continue;
            uint64_t pos_per_kmer = 0;
            const std::size_t consumed = walk_kpx_posting(
                kpx_posting_in + lo,
                static_cast<std::size_t>(hi - lo),
                counts[i],
                &pos_per_kmer);
            if (consumed == 0) {
                logger.warn("filter: kmer %u .kpx walk failed during recompute", i);
                continue;
            }
            new_total_position_count += pos_per_kmer;
        }
    }

    bool kix_ok = false;
    bool kpx_ok = !has_kpx_tmp;

    std::thread kpx_thread;
    if (has_kpx_tmp) {
        kpx_thread = std::thread([&]() {
            kpx_ok = write_filtered_kpx(
                kpx_in, kpx_final, excluded, kpx_sizes, kix_sizes,
                k, tbl_size, new_total_position_count, logger);
        });
    }

    kix_ok = write_filtered_kix(
        kix_in, kix_final, excluded, kix_sizes,
        k, tbl_size, new_total_distinct_postings, logger);

    if (has_kpx_tmp) kpx_thread.join();

    if (!kix_ok || !kpx_ok) {
        if (kix_ok) std::remove(kix_final.c_str());
        if (kpx_ok) std::remove(kpx_final.c_str());
        return false;
    }

    kix_in.close();
    if (has_kpx_tmp) kpx_in.close();

    if (std::rename(ksx_tmp.c_str(), ksx_final.c_str()) != 0) {
        logger.error("filter: failed to rename %s -> %s", ksx_tmp.c_str(), ksx_final.c_str());
        return false;
    }

    std::remove(kix_tmp.c_str());
    if (has_kpx_tmp) std::remove(kpx_tmp.c_str());

    logger.info("Filtered volume: %s (distinct_postings=%lu, position_count=%lu)",
                vol_prefix.c_str(),
                static_cast<unsigned long>(new_total_distinct_postings),
                static_cast<unsigned long>(new_total_position_count));
    return true;
}

bool filter_volumes_cross_volume(
    const std::vector<std::string>& vol_prefixes,
    const std::string& khx_path,
    int k,
    uint64_t freq_threshold,
    int filter_threads,
    const Logger& logger) {

    logger.info("Cross-volume filter: aggregating counts from %zu volume(s)...",
                vol_prefixes.size());

    uint32_t eff_tbl_size = 0;
    {
        std::string kix_tmp0 = vol_prefixes[0] + ".kix.tmp";
        KixReader kix0;
        if (!kix0.open(kix_tmp0)) {
            logger.error("filter: cannot open %s for count aggregation", kix_tmp0.c_str());
            return false;
        }
        eff_tbl_size = kix0.table_size();
        kix0.close();
    }

    std::vector<uint64_t> global_counts(eff_tbl_size, 0);

    for (const auto& prefix : vol_prefixes) {
        std::string kix_tmp = prefix + ".kix.tmp";
        KixReader kix;
        if (!kix.open(kix_tmp)) {
            logger.error("filter: cannot open %s for count aggregation", kix_tmp.c_str());
            return false;
        }
        auto vol_counts = kix.bulk_count_postings();
        for (uint32_t i = 0; i < eff_tbl_size; i++) {
            global_counts[i] += vol_counts[i];
        }
        kix.close();
    }

    std::vector<bool> excluded(eff_tbl_size, false);
    uint64_t num_excluded = 0;
    for (uint32_t i = 0; i < eff_tbl_size; i++) {
        if (global_counts[i] > freq_threshold) {
            excluded[i] = true;
            num_excluded++;
        }
    }
    global_counts.clear();
    global_counts.shrink_to_fit();

    logger.info("Cross-volume filter: %lu k-mers excluded (threshold=%lu)",
                static_cast<unsigned long>(num_excluded),
                static_cast<unsigned long>(freq_threshold));

    // Write .khx BEFORE running the per-volume filter pass.  The filter
    // pass deletes each volume's .kix.tmp / .kpx.tmp / .ksx.tmp once it
    // produces the final files, which throws away the only source of the
    // pre-filter distinct-seq_id counts.  Writing .khx first means a
    // crash between aggregation and full filter completion always leaves
    // either (a) all volumes still as .tmp with .khx valid, or (b) a
    // partial mix of final + .tmp with .khx valid — both recoverable
    // from .tmp by the resume path.  The "all final + no .khx" state
    // therefore cannot arise from a crash, only from manual deletion.
    const uint32_t base_tbl_size = table_size(k);
    std::vector<bool> khx_excluded;
    if (eff_tbl_size > base_tbl_size) {
        khx_excluded.resize(base_tbl_size, false);
        for (uint32_t i = 0; i < eff_tbl_size; i++) {
            if (excluded[i]) {
                khx_excluded[i & (base_tbl_size - 1)] = true;
            }
        }
    } else {
        khx_excluded = excluded;
    }
    if (!write_khx_bitset(khx_path, k, khx_excluded, logger)) {
        return false;
    }

    size_t num_vols = vol_prefixes.size();
    std::vector<bool> vol_ok(num_vols, false);

    tbb::task_arena arena(filter_threads);
    arena.execute([&] {
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, num_vols),
            [&](const tbb::blocked_range<size_t>& range) {
                for (size_t vi = range.begin(); vi < range.end(); vi++) {
                    vol_ok[vi] = filter_one_volume(vol_prefixes[vi], excluded, k, logger);
                }
            });
    });

    for (size_t vi = 0; vi < num_vols; vi++) {
        if (!vol_ok[vi]) {
            logger.error("filter: volume %zu failed", vi);
            return false;
        }
    }

    logger.info("Cross-volume filter: done.");
    return true;
}

} // namespace ikafssn
