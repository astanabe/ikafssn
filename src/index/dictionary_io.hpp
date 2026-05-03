#pragma once

// Phase 7a + 7e — Shared helpers for emitting the .kix / .kpx
// dictionaries as Elias-Fano blobs.  Used by KixWriter / KpxWriter,
// index_builder.cpp finalize, and index_filter.cpp.

#include <cstdint>
#include <cstdio>

namespace ikafssn {

// Encode the .kix dictionary (offsets[0..table_size+1] sentinel-
// terminated, monotonically non-decreasing) as Elias-Fano and write
// the resulting blob to `fp`.  Returns false on I/O error.
//
//   offsets             length = table_size + 1
//   table_size          = 4^k (entry count, sentinel at offsets[table_size])
//   posting_file_size   exclusive upper bound on offset values (= the
//                       sentinel value); used as U_raw for the EF
//                       universe so the encoder does not have to scan.
bool write_kix_dictionary_ef(std::FILE* fp,
                             const std::uint64_t* offsets,
                             std::uint32_t table_size,
                             std::uint64_t posting_file_size);

// Encode the .kpx pos_offsets dictionary as Elias-Fano and write the
// blob to `fp`.  Like .kix but without a sentinel — the .kpx posting
// lists are self-delimiting and callers use posting_file_size as a
// loose upper bound when fetching the last k-mer's posting list.
//
//   pos_offsets         length = table_size
//   table_size          = 4^k
//   posting_file_size   exclusive upper bound on pos_offset values; used
//                       as U_raw for the EF universe.
bool write_kpx_dictionary_ef(std::FILE* fp,
                             const std::uint64_t* pos_offsets,
                             std::uint32_t table_size,
                             std::uint64_t posting_file_size);

} // namespace ikafssn
