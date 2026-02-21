#include "index/index_filter.hpp"
#include "index/kix_reader.hpp"
#include "index/kpx_reader.hpp"
#include "index/kix_format.hpp"
#include "index/kpx_format.hpp"
#include "index/khx_writer.hpp"
#include "core/config.hpp"
#include "core/types.hpp"
#include "util/logger.hpp"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <thread>
#include <vector>

#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

namespace ikafssn {

// Compute the byte size of each k-mer's posting data from offsets and counts.
// For k-mers with count==0, size is 0.
// For the last non-zero k-mer, its data extends to total_data_size.
static std::vector<uint64_t> compute_posting_sizes(
    const uint64_t* offsets, const uint32_t* counts,
    uint64_t tbl_size, uint64_t total_data_size) {

    std::vector<uint64_t> sizes(tbl_size, 0);

    // Find all k-mers with postings, in order
    uint64_t prev_kmer = UINT64_MAX;
    uint64_t prev_offset = 0;

    for (uint64_t i = 0; i < tbl_size; i++) {
        if (counts[i] > 0) {
            if (prev_kmer != UINT64_MAX) {
                sizes[prev_kmer] = offsets[i] - prev_offset;
            }
            prev_kmer = i;
            prev_offset = offsets[i];
        }
    }
    if (prev_kmer != UINT64_MAX) {
        sizes[prev_kmer] = total_data_size - prev_offset;
    }

    return sizes;
}

// Write filtered .kix from .kix.tmp (excluded k-mers removed).
// Returns true on success.
static bool write_filtered_kix(
    const KixReader& kix_in,
    const std::string& kix_final,
    const std::vector<bool>& excluded,
    const std::vector<uint64_t>& kix_sizes,
    int k,
    uint64_t new_total_postings,
    const Logger& logger) {

    const uint64_t tbl_size = table_size(k);
    const uint64_t* kix_offsets_in = kix_in.offsets();
    const uint32_t* counts_in = kix_in.counts();
    const uint8_t* kix_posting_in = kix_in.posting_data();

    FILE* kix_fp = std::fopen(kix_final.c_str(), "wb");
    if (!kix_fp) {
        logger.error("filter: cannot open %s for writing", kix_final.c_str());
        return false;
    }

    // Write header placeholder
    KixHeader kix_hdr{};
    std::fwrite(&kix_hdr, sizeof(kix_hdr), 1, kix_fp);

    // Build new offsets and counts, write posting data
    std::vector<uint64_t> new_kix_offsets(tbl_size, 0);
    std::vector<uint32_t> new_counts(tbl_size, 0);

    // Reserve space for offsets table
    std::fwrite(new_kix_offsets.data(), sizeof(uint64_t), tbl_size, kix_fp);
    // Reserve space for counts table
    std::fwrite(new_counts.data(), sizeof(uint32_t), tbl_size, kix_fp);

    uint64_t kix_data_pos = 0;
    for (uint64_t i = 0; i < tbl_size; i++) {
        if (counts_in[i] > 0 && !excluded[i]) {
            new_kix_offsets[i] = kix_data_pos;
            new_counts[i] = counts_in[i];
            std::fwrite(kix_posting_in + kix_offsets_in[i], 1, kix_sizes[i], kix_fp);
            kix_data_pos += kix_sizes[i];
        }
    }

    // Write header
    std::memcpy(kix_hdr.magic, KIX_MAGIC, 4);
    kix_hdr.format_version = KIX_FORMAT_VERSION;
    kix_hdr.k = static_cast<uint8_t>(k);
    kix_hdr.kmer_type = kmer_type_for_k(k);
    kix_hdr.num_sequences = kix_in.num_sequences();
    kix_hdr.total_postings = new_total_postings;
    kix_hdr.flags = kix_in.header().flags;
    kix_hdr.volume_index = kix_in.header().volume_index;
    kix_hdr.total_volumes = kix_in.header().total_volumes;
    kix_hdr.db_name_len = kix_in.header().db_name_len;
    std::memcpy(kix_hdr.db_name, kix_in.header().db_name, 32);

    std::fseek(kix_fp, 0, SEEK_SET);
    std::fwrite(&kix_hdr, sizeof(kix_hdr), 1, kix_fp);
    std::fwrite(new_kix_offsets.data(), sizeof(uint64_t), tbl_size, kix_fp);
    std::fwrite(new_counts.data(), sizeof(uint32_t), tbl_size, kix_fp);

    std::fclose(kix_fp);
    return true;
}

// Write filtered .kpx from .kpx.tmp (excluded k-mers removed).
// Returns true on success.
static bool write_filtered_kpx(
    const KpxReader& kpx_in,
    const uint32_t* counts_in,
    const std::string& kpx_final,
    const std::vector<bool>& excluded,
    const std::vector<uint64_t>& kpx_sizes,
    int k,
    uint64_t new_total_postings,
    const Logger& logger) {

    const uint64_t tbl_size = table_size(k);
    const uint64_t* kpx_offsets_in = kpx_in.pos_offsets();
    const uint8_t* kpx_posting_in = kpx_in.posting_data();

    FILE* kpx_fp = std::fopen(kpx_final.c_str(), "wb");
    if (!kpx_fp) {
        logger.error("filter: cannot open %s for writing", kpx_final.c_str());
        return false;
    }

    // Write header placeholder
    KpxHeader kpx_hdr{};
    std::fwrite(&kpx_hdr, sizeof(kpx_hdr), 1, kpx_fp);

    std::vector<uint64_t> new_kpx_offsets(tbl_size, 0);
    std::fwrite(new_kpx_offsets.data(), sizeof(uint64_t), tbl_size, kpx_fp);

    uint64_t kpx_data_pos = 0;
    for (uint64_t i = 0; i < tbl_size; i++) {
        if (counts_in[i] > 0 && !excluded[i]) {
            new_kpx_offsets[i] = kpx_data_pos;
            std::fwrite(kpx_posting_in + kpx_offsets_in[i], 1, kpx_sizes[i], kpx_fp);
            kpx_data_pos += kpx_sizes[i];
        }
    }

    // Write header
    std::memcpy(kpx_hdr.magic, KPX_MAGIC, 4);
    kpx_hdr.format_version = KPX_FORMAT_VERSION;
    kpx_hdr.k = static_cast<uint8_t>(k);
    kpx_hdr.total_postings = new_total_postings;

    std::fseek(kpx_fp, 0, SEEK_SET);
    std::fwrite(&kpx_hdr, sizeof(kpx_hdr), 1, kpx_fp);
    std::fwrite(new_kpx_offsets.data(), sizeof(uint64_t), tbl_size, kpx_fp);

    std::fclose(kpx_fp);
    return true;
}

// Filter a single volume's .kix.tmp/.kpx.tmp -> .kix/.kpx (in parallel).
// Also renames .ksx.tmp -> .ksx and removes .tmp files.
static bool filter_one_volume(
    const std::string& vol_prefix,
    const std::vector<bool>& excluded,
    int k,
    const Logger& logger) {

    const uint64_t tbl_size = table_size(k);

    std::string kix_tmp = vol_prefix + ".kix.tmp";
    std::string kpx_tmp = vol_prefix + ".kpx.tmp";
    std::string ksx_tmp = vol_prefix + ".ksx.tmp";
    std::string kix_final = vol_prefix + ".kix";
    std::string kpx_final = vol_prefix + ".kpx";
    std::string ksx_final = vol_prefix + ".ksx";

    // Open .kix.tmp and .kpx.tmp via readers (mmap)
    KixReader kix_in;
    if (!kix_in.open(kix_tmp)) {
        logger.error("filter: cannot open %s", kix_tmp.c_str());
        return false;
    }
    KpxReader kpx_in;
    if (!kpx_in.open(kpx_tmp)) {
        logger.error("filter: cannot open %s", kpx_tmp.c_str());
        return false;
    }

    const uint32_t* counts_in = kix_in.counts();

    // Compute posting byte sizes
    auto kix_sizes = compute_posting_sizes(
        kix_in.offsets(), counts_in, tbl_size, kix_in.posting_data_size());
    auto kpx_sizes = compute_posting_sizes(
        kpx_in.pos_offsets(), counts_in, tbl_size, kpx_in.posting_data_size());

    // Compute new totals
    uint64_t new_total_postings = 0;
    for (uint64_t i = 0; i < tbl_size; i++) {
        if (!excluded[i]) {
            new_total_postings += counts_in[i];
        }
    }

    // Write .kix and .kpx in parallel (two threads)
    bool kix_ok = false;
    bool kpx_ok = false;

    std::thread kpx_thread([&]() {
        kpx_ok = write_filtered_kpx(
            kpx_in, counts_in, kpx_final, excluded, kpx_sizes,
            k, new_total_postings, logger);
    });

    kix_ok = write_filtered_kix(
        kix_in, kix_final, excluded, kix_sizes,
        k, new_total_postings, logger);

    kpx_thread.join();

    if (!kix_ok || !kpx_ok) {
        if (kix_ok) std::remove(kix_final.c_str());
        if (kpx_ok) std::remove(kpx_final.c_str());
        return false;
    }

    // Close mmapped readers before removing .tmp files
    kix_in.close();
    kpx_in.close();

    // Rename .ksx.tmp -> .ksx
    if (std::rename(ksx_tmp.c_str(), ksx_final.c_str()) != 0) {
        logger.error("filter: failed to rename %s -> %s", ksx_tmp.c_str(), ksx_final.c_str());
        return false;
    }

    // Remove .tmp files
    std::remove(kix_tmp.c_str());
    std::remove(kpx_tmp.c_str());

    logger.info("Filtered volume: %s (total_postings: %lu)",
                vol_prefix.c_str(), static_cast<unsigned long>(new_total_postings));
    return true;
}

bool filter_volumes_cross_volume(
    const std::vector<std::string>& vol_prefixes,
    const std::string& khx_path,
    int k,
    uint64_t freq_threshold,
    const Logger& logger) {

    const uint64_t tbl_size = table_size(k);

    // Step 1: Aggregate counts across all volumes
    logger.info("Cross-volume filter: aggregating counts from %zu volume(s)...",
                vol_prefixes.size());

    std::vector<uint64_t> global_counts(tbl_size, 0);

    for (const auto& prefix : vol_prefixes) {
        std::string kix_tmp = prefix + ".kix.tmp";
        KixReader kix;
        if (!kix.open(kix_tmp)) {
            logger.error("filter: cannot open %s for count aggregation", kix_tmp.c_str());
            return false;
        }
        const uint32_t* cts = kix.counts();
        for (uint64_t i = 0; i < tbl_size; i++) {
            global_counts[i] += cts[i];
        }
        kix.close();
    }

    // Step 2: Determine excluded k-mers
    std::vector<bool> excluded(tbl_size, false);
    uint64_t num_excluded = 0;
    for (uint64_t i = 0; i < tbl_size; i++) {
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

    // Step 3: Filter each volume in parallel
    size_t num_vols = vol_prefixes.size();
    std::vector<bool> vol_ok(num_vols, false);

    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, num_vols),
        [&](const tbb::blocked_range<size_t>& range) {
            for (size_t vi = range.begin(); vi < range.end(); vi++) {
                vol_ok[vi] = filter_one_volume(vol_prefixes[vi], excluded, k, logger);
            }
        });

    for (size_t vi = 0; vi < num_vols; vi++) {
        if (!vol_ok[vi]) {
            logger.error("filter: volume %zu failed", vi);
            return false;
        }
    }

    // Step 4: Write shared .khx
    if (!write_khx_bitset(khx_path, k, excluded, logger)) {
        return false;
    }

    logger.info("Cross-volume filter: done.");
    return true;
}

} // namespace ikafssn
