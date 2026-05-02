#pragma once

// Phase 7a — Shared helper for emitting the .kix dictionary as an
// Elias-Fano blob.  Used by KixWriter, index_builder.cpp finalize, and
// index_filter.cpp::write_filtered_kix to keep dictionary emission in
// one place.
//
// The .kpx pos_offsets dictionary keeps the v8 raw u32/u64 layout in
// 7a; Phase 7e moves it to EF as well.

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

} // namespace ikafssn
