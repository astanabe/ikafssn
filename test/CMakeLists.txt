function(add_ikafssn_test TEST_NAME TEST_SOURCE)
    add_executable(${TEST_NAME} ${TEST_SOURCE})
    target_link_libraries(${TEST_NAME} PRIVATE ikafssn_index ikafssn_io ikafssn_core)
    target_include_directories(${TEST_NAME} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endfunction()

# Register SIMD-aware variants of a test. For each tier (x86_64 + aarch64),
# a separate test invocation is created with IKAFSSN_FORCE_SIMD pinning the
# kernel selection and IKAFSSN_TEST_REQUIRE_TIER ensuring CPUs that do not
# actually support that tier exit(77) (CTest skip).
function(add_ikafssn_simd_test TEST_NAME)
    # x86_64 4-tier ladder (SSE4.2 floor, no standalone VBMI).
    # The Scalar tier is intentionally omitted: init_simd_dispatch()
    # rejects auto_cap == Scalar with exit(2) so a `_force_scalar`
    # variant would always fail rather than skip cleanly.
    set(_simd_tiers
        sse42 avx2 avx512bw avx512vbmi2
        neon sve sve2
    )
    foreach(tier IN LISTS _simd_tiers)
        add_test(NAME ${TEST_NAME}_force_${tier} COMMAND ${TEST_NAME})
        set_tests_properties(${TEST_NAME}_force_${tier} PROPERTIES
            ENVIRONMENT "IKAFSSN_FORCE_SIMD=${tier};IKAFSSN_TEST_REQUIRE_TIER=${tier}"
            SKIP_RETURN_CODE 77
            LABELS "simd_variant"
        )
    endforeach()
endfunction()

add_ikafssn_test(test_varint test_varint.cpp)

# Compressed I/O round-trip + concatenation + magic / extension detection.
# ikafssn_io already PRIVATE-links zlib/bz2/lzma/zstd, so the test
# binary inherits them via the static archive.
add_ikafssn_test(test_compressed_stream test_compressed_stream.cpp)

add_ikafssn_test(test_kmer_encoding test_kmer_encoding.cpp)
target_link_libraries(test_kmer_encoding PRIVATE ikafssn_util)
target_compile_definitions(test_kmer_encoding PRIVATE IKAFSSN_TESTING)
add_ikafssn_simd_test(test_kmer_encoding)

# FASTA toupper SIMD test
add_ikafssn_test(test_fasta_toupper test_fasta_toupper.cpp)
target_link_libraries(test_fasta_toupper PRIVATE ikafssn_util)
target_compile_definitions(test_fasta_toupper PRIVATE IKAFSSN_TESTING)
add_ikafssn_simd_test(test_fasta_toupper)

# ncbi2na unpack SIMD test
add_ikafssn_test(test_ncbi2na_unpack test_ncbi2na_unpack.cpp)
target_link_libraries(test_ncbi2na_unpack PRIVATE ikafssn_util)
target_compile_definitions(test_ncbi2na_unpack PRIVATE IKAFSSN_TESTING)
add_ikafssn_simd_test(test_ncbi2na_unpack)

# SIMD dispatch unit test (no external dependencies).
add_executable(test_simd_dispatch test_simd_dispatch.cpp)
target_link_libraries(test_simd_dispatch PRIVATE ikafssn_util ikafssn_core)
target_include_directories(test_simd_dispatch PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_compile_definitions(test_simd_dispatch PRIVATE IKAFSSN_TESTING)
add_test(NAME test_simd_dispatch COMMAND test_simd_dispatch)
add_ikafssn_simd_test(test_simd_dispatch)
add_ikafssn_test(test_spaced_seed test_spaced_seed.cpp)
target_link_libraries(test_spaced_seed PRIVATE ikafssn_util)
target_compile_definitions(test_spaced_seed PRIVATE IKAFSSN_TESTING)
add_ikafssn_simd_test(test_spaced_seed)

# extract_for_mask_batch SIMD kernel test
add_ikafssn_test(test_extract_for_mask_batch test_extract_for_mask_batch.cpp)
target_link_libraries(test_extract_for_mask_batch PRIVATE ikafssn_util)
target_compile_definitions(test_extract_for_mask_batch PRIVATE IKAFSSN_TESTING)
add_ikafssn_simd_test(test_extract_for_mask_batch)
add_ikafssn_test(test_ksx_io test_ksx_io.cpp)
add_ikafssn_test(test_kix_io test_kix_io.cpp)
add_ikafssn_test(test_kpx_io test_kpx_io.cpp)
add_ikafssn_test(test_khx_io test_khx_io.cpp)
target_link_libraries(test_khx_io PRIVATE ikafssn_util)

# fragment_splitter unit test (no BLAST DB dependency)
add_ikafssn_test(test_fragment_splitter test_fragment_splitter.cpp)

# result_dedup unit test (uses KsxWriter/Reader for the Stage 2
# parent-relative key resolution; no BLAST DB dependency).
add_executable(test_result_dedup test_result_dedup.cpp)
target_link_libraries(test_result_dedup PRIVATE
    ikafssn_search ikafssn_index ikafssn_io ikafssn_core
)
target_include_directories(test_result_dedup PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
)
add_test(NAME test_result_dedup COMMAND test_result_dedup)

# SIMD dedup kernel unit test
add_ikafssn_test(test_seq_id_dedup test_seq_id_dedup.cpp)
target_link_libraries(test_seq_id_dedup PRIVATE ikafssn_util)
target_compile_definitions(test_seq_id_dedup PRIVATE IKAFSSN_TESTING)
add_ikafssn_simd_test(test_seq_id_dedup)

# Elias-Fano dictionary codec round-trip test, exercised via every
# per-ISA tier object (sse42 / avx2 / avx512bw / avx512vbmi2 / neon)
# by add_ikafssn_simd_test.
add_ikafssn_test(test_ef_codec test_ef_codec.cpp)
target_link_libraries(test_ef_codec PRIVATE ikafssn_util)
target_compile_definitions(test_ef_codec PRIVATE IKAFSSN_TESTING)
add_ikafssn_simd_test(test_ef_codec)

# Confirm the EF dispatcher selects the tier requested via
# IKAFSSN_FORCE_SIMD (mirrors test_simd_dispatch's approach).
add_ikafssn_test(test_ef_simd_dispatch test_ef_simd_dispatch.cpp)
target_link_libraries(test_ef_simd_dispatch PRIVATE ikafssn_util)
target_compile_definitions(test_ef_simd_dispatch PRIVATE IKAFSSN_TESTING)
add_ikafssn_simd_test(test_ef_simd_dispatch)

# Synthesize older-format index files and verify that current readers
# reject them with the expected error message.
add_ikafssn_test(test_v8_rejection test_v8_rejection.cpp)
target_link_libraries(test_v8_rejection PRIVATE ikafssn_util)

# Structural validation of (.kix, .kpx) volume pairs + bit-flip
# injection regression.  Exercises validate_volume() against both a
# clean build and a corrupted kind map.
add_ikafssn_test(test_index_filter_validate test_index_filter_validate.cpp)
target_link_libraries(test_index_filter_validate PRIVATE ikafssn_util)

# ikafssnindex resume / per-volume validation regression.  Spawns the
# ikafssnindex CLI binary against the SSU fixture.
add_executable(test_resume_validation test_resume_validation.cpp)
target_link_libraries(test_resume_validation PRIVATE
    ikafssn_index ikafssn_io ikafssn_util ikafssn_core
)
target_include_directories(test_resume_validation PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_SOURCE_DIR}/src"
)
target_compile_definitions(test_resume_validation PRIVATE
    SOURCE_DIR="${CMAKE_SOURCE_DIR}"
)
add_test(NAME test_resume_validation COMMAND test_resume_validation)

# Per-tier popcount_kinds + horizontal_sum_u8 unit test.  Run across
# every ISA tier via add_ikafssn_simd_test so AVX2 / AVX-512
# auto-vectorisation paths get the same coverage as SSE4.2 / NEON.
add_ikafssn_test(test_kind_map_popcount test_kind_map_popcount.cpp)
target_link_libraries(test_kind_map_popcount PRIVATE ikafssn_util)
target_compile_definitions(test_kind_map_popcount PRIVATE IKAFSSN_TESTING)
add_ikafssn_simd_test(test_kind_map_popcount)

# Builder integration test (requires NCBI Toolkit + BLAST DB)
add_executable(test_builder test_builder.cpp)
target_link_libraries(test_builder PRIVATE
    ikafssn_builder ikafssn_index ikafssn_blastdb ikafssn_io ikafssn_util ikafssn_core
)
target_include_directories(test_builder PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${NCBI_TOOLKIT_INCLUDE}"
)
target_compile_definitions(test_builder PRIVATE
    SOURCE_DIR="${CMAKE_SOURCE_DIR}"
    IKAFSSN_TESTING
)
add_test(NAME test_builder COMMAND test_builder)
add_ikafssn_simd_test(test_builder)

# parallel_sort_dispatch unit test (no BLAST DB dependency)
add_executable(test_parallel_sort test_parallel_sort.cpp)
target_link_libraries(test_parallel_sort PRIVATE
    ikafssn_index ikafssn_io ikafssn_util ikafssn_core
)
target_include_directories(test_parallel_sort PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
)
target_compile_definitions(test_parallel_sort PRIVATE IKAFSSN_TESTING)
add_test(NAME test_parallel_sort COMMAND test_parallel_sort)
add_ikafssn_simd_test(test_parallel_sort)

# PackedKmerScanner + AmbiguityParser test (requires NCBI Toolkit + BLAST DB)
add_executable(test_packed_kmer_scanner test_packed_kmer_scanner.cpp)
target_link_libraries(test_packed_kmer_scanner PRIVATE
    ikafssn_builder ikafssn_index ikafssn_blastdb ikafssn_io ikafssn_util ikafssn_core
)
target_include_directories(test_packed_kmer_scanner PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${NCBI_TOOLKIT_INCLUDE}"
)
target_compile_definitions(test_packed_kmer_scanner PRIVATE
    SOURCE_DIR="${CMAKE_SOURCE_DIR}"
    IKAFSSN_TESTING
)
add_test(NAME test_packed_kmer_scanner COMMAND test_packed_kmer_scanner)
add_ikafssn_simd_test(test_packed_kmer_scanner)

# --- Search-pipeline tests ---

# Diagonal filter test (no external dependencies)
add_executable(test_diagonal_filter test_diagonal_filter.cpp)
target_link_libraries(test_diagonal_filter PRIVATE
    ikafssn_search ikafssn_index ikafssn_io ikafssn_core
)
target_include_directories(test_diagonal_filter PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
add_test(NAME test_diagonal_filter COMMAND test_diagonal_filter)

# Chaining DP test (no external dependencies)
add_executable(test_chaining test_chaining.cpp)
target_link_libraries(test_chaining PRIVATE
    ikafssn_search ikafssn_index ikafssn_io ikafssn_core
)
target_include_directories(test_chaining PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
add_test(NAME test_chaining COMMAND test_chaining)

# Seqidlist reader test (no external dependencies)
add_executable(test_seqidlist_reader test_seqidlist_reader.cpp)
target_link_libraries(test_seqidlist_reader PRIVATE
    ikafssn_io ikafssn_core
)
target_include_directories(test_seqidlist_reader PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
add_test(NAME test_seqidlist_reader COMMAND test_seqidlist_reader)

# OID filter test (no external dependencies, uses KsxWriter/Reader)
add_executable(test_oid_filter test_oid_filter.cpp)
target_link_libraries(test_oid_filter PRIVATE
    ikafssn_search ikafssn_index ikafssn_io ikafssn_core
)
target_include_directories(test_oid_filter PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
add_test(NAME test_oid_filter COMMAND test_oid_filter)

# Stage 1 filter test (requires NCBI Toolkit for index building)
add_executable(test_stage1 test_stage1.cpp)
target_link_libraries(test_stage1 PRIVATE
    ikafssn_search ikafssn_builder ikafssn_index ikafssn_blastdb ikafssn_io ikafssn_util ikafssn_core
)
target_include_directories(test_stage1 PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${NCBI_TOOLKIT_INCLUDE}"
)
target_compile_definitions(test_stage1 PRIVATE
    SOURCE_DIR="${CMAKE_SOURCE_DIR}"
    IKAFSSN_TESTING
)
add_test(NAME test_stage1 COMMAND test_stage1)
add_ikafssn_simd_test(test_stage1)

# Stage 1 SIMD kernel unit test (no BLAST DB dependency)
add_executable(test_stage1_simd_kernel test_stage1_simd_kernel.cpp)
target_link_libraries(test_stage1_simd_kernel PRIVATE
    ikafssn_search ikafssn_index ikafssn_io ikafssn_util ikafssn_core
)
target_include_directories(test_stage1_simd_kernel PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
)
target_compile_definitions(test_stage1_simd_kernel PRIVATE IKAFSSN_TESTING)
add_test(NAME test_stage1_simd_kernel COMMAND test_stage1_simd_kernel)
add_ikafssn_simd_test(test_stage1_simd_kernel)

# Split-index end-to-end integration test (requires NCBI Toolkit +
# BLAST DB + the long-parent fixture from setup_ssu_testdata.sh).
add_executable(test_split_index_integration test_split_index_integration.cpp)
target_link_libraries(test_split_index_integration PRIVATE
    ikafssn_search ikafssn_builder ikafssn_index ikafssn_blastdb ikafssn_io ikafssn_util ikafssn_core
)
target_include_directories(test_split_index_integration PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${NCBI_TOOLKIT_INCLUDE}"
    "${CMAKE_SOURCE_DIR}/src"
)
target_compile_definitions(test_split_index_integration PRIVATE
    SOURCE_DIR="${CMAKE_SOURCE_DIR}"
)
add_test(NAME test_split_index_integration COMMAND test_split_index_integration)

# Search integration test (requires NCBI Toolkit + BLAST DB)
add_executable(test_search_integration test_search_integration.cpp)
target_link_libraries(test_search_integration PRIVATE
    ikafssn_search ikafssn_builder ikafssn_index ikafssn_blastdb ikafssn_io ikafssn_util ikafssn_core
)
target_include_directories(test_search_integration PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${NCBI_TOOLKIT_INCLUDE}"
)
target_compile_definitions(test_search_integration PRIVATE
    SOURCE_DIR="${CMAKE_SOURCE_DIR}"
)
add_test(NAME test_search_integration COMMAND test_search_integration)

# --- Multi-volume tests ---

# Multi-volume and parallel search test (requires NCBI Toolkit + TBB)
add_executable(test_multivolume test_multivolume.cpp)
target_link_libraries(test_multivolume PRIVATE
    ikafssn_search ikafssn_builder ikafssn_index ikafssn_blastdb ikafssn_io ikafssn_util ikafssn_core
    TBB::tbb
)
target_include_directories(test_multivolume PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${NCBI_TOOLKIT_INCLUDE}"
)
target_compile_definitions(test_multivolume PRIVATE
    SOURCE_DIR="${CMAKE_SOURCE_DIR}"
)
add_test(NAME test_multivolume COMMAND test_multivolume)

# --- Server / protocol tests ---

# Protocol serialize/deserialize test (no external dependencies)
add_executable(test_protocol test_protocol.cpp)
target_link_libraries(test_protocol PRIVATE
    ikafssn_protocol ikafssn_util ikafssn_core
)
target_include_directories(test_protocol PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
add_test(NAME test_protocol COMMAND test_protocol)

# Inter-request posting-budget pool unit test (no server, no BLAST DB)
find_package(Threads REQUIRED)
add_executable(test_budget_pool test_budget_pool.cpp
    "${CMAKE_SOURCE_DIR}/src/ikafssnserver/budget_pool.cpp"
)
target_include_directories(test_budget_pool PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_SOURCE_DIR}/src"
)
target_link_libraries(test_budget_pool PRIVATE Threads::Threads)
add_test(NAME test_budget_pool COMMAND test_budget_pool)

# Server-client integration test (requires NCBI Toolkit + TBB)
add_executable(test_server_client test_server_client.cpp)
target_link_libraries(test_server_client PRIVATE
    ikafssn_server_lib ikafssn_client_lib
    ikafssn_search ikafssn_builder ikafssn_index ikafssn_blastdb ikafssn_io ikafssn_util ikafssn_core
    TBB::tbb
)
target_include_directories(test_server_client PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${NCBI_TOOLKIT_INCLUDE}"
    "${CMAKE_SOURCE_DIR}/src"
)
target_compile_definitions(test_server_client PRIVATE
    SOURCE_DIR="${CMAKE_SOURCE_DIR}"
)
add_test(NAME test_server_client COMMAND test_server_client)

# --- HTTPD tests ---

# HTTPD backend client integration test (requires NCBI Toolkit + TBB + Drogon build)
if(BUILD_HTTPD)
    add_executable(test_httpd test_httpd.cpp)
    target_link_libraries(test_httpd PRIVATE
        ikafssn_httpd_backend ikafssn_server_lib
        ikafssn_search ikafssn_builder ikafssn_index ikafssn_blastdb ikafssn_io ikafssn_util ikafssn_core
        TBB::tbb
    )
    target_include_directories(test_httpd PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}"
        "${NCBI_TOOLKIT_INCLUDE}"
        "${CMAKE_SOURCE_DIR}/src"
    )
    target_compile_definitions(test_httpd PRIVATE
        SOURCE_DIR="${CMAKE_SOURCE_DIR}"
    )
    add_test(NAME test_httpd COMMAND test_httpd)

    # Full stack integration test: ikafssnclient (HTTP) -> ikafssnhttpd -> ikafssnserver
    if(ENABLE_CLIENT_HTTP)
        add_executable(test_httpd_client test_httpd_client.cpp)
        target_link_libraries(test_httpd_client PRIVATE
            ikafssn_client_lib ikafssn_httpd_controller ikafssn_server_lib
            ikafssn_search ikafssn_builder ikafssn_index ikafssn_blastdb ikafssn_io ikafssn_util ikafssn_core
            TBB::tbb
        )
        target_include_directories(test_httpd_client PRIVATE
            "${CMAKE_CURRENT_SOURCE_DIR}"
            "${NCBI_TOOLKIT_INCLUDE}"
            "${CMAKE_SOURCE_DIR}/src"
        )
        target_compile_definitions(test_httpd_client PRIVATE
            SOURCE_DIR="${CMAKE_SOURCE_DIR}"
        )
        add_test(NAME test_httpd_client COMMAND test_httpd_client)
    endif()
endif()

# --- Retrieve tests ---

# Result reader test (no external dependencies)
add_executable(test_result_reader test_result_reader.cpp)
target_link_libraries(test_result_reader PRIVATE
    ikafssn_io ikafssn_core
)
target_include_directories(test_result_reader PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
add_test(NAME test_result_reader COMMAND test_result_reader)

# efetch retriever unit test (no external dependencies, no HTTP calls)
add_executable(test_efetch_retriever test_efetch_retriever.cpp
    "${CMAKE_SOURCE_DIR}/src/ikafssnretrieve/efetch_retriever.cpp"
)
target_link_libraries(test_efetch_retriever PRIVATE
    ikafssn_io ikafssn_core
)
target_include_directories(test_efetch_retriever PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_SOURCE_DIR}/src"
)
add_test(NAME test_efetch_retriever COMMAND test_efetch_retriever)

# Retrieve pipeline integration test (requires NCBI Toolkit)
add_executable(test_retrieve_pipeline test_retrieve_pipeline.cpp
    "${CMAKE_SOURCE_DIR}/src/ikafssnretrieve/local_retriever.cpp"
)
target_link_libraries(test_retrieve_pipeline PRIVATE
    ikafssn_search ikafssn_builder ikafssn_index ikafssn_blastdb ikafssn_io ikafssn_util ikafssn_core
)
target_include_directories(test_retrieve_pipeline PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${NCBI_TOOLKIT_INCLUDE}"
    "${CMAKE_SOURCE_DIR}/src"
)
target_compile_definitions(test_retrieve_pipeline PRIVATE
    SOURCE_DIR="${CMAKE_SOURCE_DIR}"
)
add_test(NAME test_retrieve_pipeline COMMAND test_retrieve_pipeline)

# SSU_eukaryote_rRNA integration test (requires network + NCBI Toolkit; optionally libcurl)
add_executable(test_retrieve_ssu test_retrieve_ssu.cpp
    "${CMAKE_SOURCE_DIR}/src/ikafssnretrieve/local_retriever.cpp"
    "${CMAKE_SOURCE_DIR}/src/ikafssnretrieve/efetch_retriever.cpp"
)
target_link_libraries(test_retrieve_ssu PRIVATE
    ikafssn_blastdb ikafssn_io ikafssn_util ikafssn_core
)
target_include_directories(test_retrieve_ssu PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${NCBI_TOOLKIT_INCLUDE}"
    "${CMAKE_SOURCE_DIR}/src"
)
if(ENABLE_REMOTE_RETRIEVE)
    target_link_libraries(test_retrieve_ssu PRIVATE CURL::libcurl)
    target_compile_definitions(test_retrieve_ssu PRIVATE IKAFSSN_ENABLE_REMOTE)
endif()
add_test(NAME test_retrieve_ssu COMMAND test_retrieve_ssu)

# --- Stage 3 tests ---

# Stage 3 alignment unit test (Parasail only, no BLAST DB)
add_executable(test_stage3_alignment test_stage3_alignment.cpp)
target_link_libraries(test_stage3_alignment PRIVATE
    ikafssn_stage3 ikafssn_search ikafssn_index ikafssn_blastdb
    ikafssn_io ikafssn_util ikafssn_core parasail TBB::tbb)
target_include_directories(test_stage3_alignment PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}" "${NCBI_TOOLKIT_INCLUDE}" "${CMAKE_SOURCE_DIR}/src")
add_test(NAME test_stage3_alignment COMMAND test_stage3_alignment)

# SAM/BAM writer test (htslib)
add_executable(test_sam_writer test_sam_writer.cpp)
target_link_libraries(test_sam_writer PRIVATE
    ikafssn_sam_writer ikafssn_io ikafssn_core hts)
target_include_directories(test_sam_writer PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_SOURCE_DIR}/src")
add_test(NAME test_sam_writer COMMAND test_sam_writer)

# Query preprocessor spec compliance (Nqkmer/Nhighfreq, skip reasons)
add_executable(test_query_preprocessor test_query_preprocessor.cpp)
target_link_libraries(test_query_preprocessor PRIVATE
    ikafssn_search ikafssn_index ikafssn_io ikafssn_protocol
    ikafssn_util ikafssn_core)
target_include_directories(test_query_preprocessor PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_SOURCE_DIR}/src")
add_test(NAME test_query_preprocessor COMMAND test_query_preprocessor)

# Skip-marker output across TSV / JSON / SAM
add_executable(test_skip_output test_skip_output.cpp)
target_link_libraries(test_skip_output PRIVATE
    ikafssn_sam_writer ikafssn_io ikafssn_protocol ikafssn_core hts)
target_include_directories(test_skip_output PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_SOURCE_DIR}/src")
add_test(NAME test_skip_output COMMAND test_skip_output)

# Stage 3 integration test (requires NCBI Toolkit + Parasail + BLAST DB)
add_executable(test_stage3_integration test_stage3_integration.cpp)
target_link_libraries(test_stage3_integration PRIVATE
    ikafssn_stage3 ikafssn_sam_writer ikafssn_search ikafssn_builder
    ikafssn_index ikafssn_blastdb ikafssn_io ikafssn_util ikafssn_core
    parasail hts TBB::tbb)
target_include_directories(test_stage3_integration PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}" "${NCBI_TOOLKIT_INCLUDE}" "${CMAKE_SOURCE_DIR}/src")
target_compile_definitions(test_stage3_integration PRIVATE
    SOURCE_DIR="${CMAKE_SOURCE_DIR}")
add_test(NAME test_stage3_integration COMMAND test_stage3_integration)

# --- Compression / job-store tests ---

# One-shot zstd compress/decompress round-trip test
if(ENABLE_CLIENT_HTTP OR BUILD_HTTPD)
    add_executable(test_zstd_io test_zstd_io.cpp)
    target_link_libraries(test_zstd_io PRIVATE
        ikafssn_zstd_oneshot ikafssn_core
    )
    target_include_directories(test_zstd_io PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_SOURCE_DIR}/src"
    )
    add_test(NAME test_zstd_io COMMAND test_zstd_io)
endif()

# JobStore round-trip + retry + TTL tests
if(BUILD_HTTPD)
    add_executable(test_job_store test_job_store.cpp)
    target_link_libraries(test_job_store PRIVATE
        ikafssn_httpd_jobstore ikafssn_core
    )
    target_include_directories(test_job_store PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_SOURCE_DIR}/src"
    )
    add_test(NAME test_job_store COMMAND test_job_store)

    # ResultStore unit test: per-job zstd file storage with shard layout.
    add_executable(test_result_store test_result_store.cpp)
    target_link_libraries(test_result_store PRIVATE
        ikafssn_httpd_jobstore ikafssn_zstd_oneshot ikafssn_core
    )
    target_include_directories(test_result_store PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_SOURCE_DIR}/src"
    )
    add_test(NAME test_result_store COMMAND test_result_store)

    # QueryStore unit test: per-job zstd request body storage,
    # passthrough vs plain write, sweep_tmp, list_job_ids, shard layout.
    add_executable(test_query_store test_query_store.cpp)
    target_link_libraries(test_query_store PRIVATE
        ikafssn_httpd_jobstore ikafssn_zstd_oneshot ikafssn_core
    )
    target_include_directories(test_query_store PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_SOURCE_DIR}/src"
    )
    add_test(NAME test_query_store COMMAND test_query_store)
endif()

# Checkpoint test (requires OpenSSL)
add_executable(test_checkpoint test_checkpoint.cpp
    "${CMAKE_SOURCE_DIR}/src/ikafssnclient/checkpoint.cpp"
)
target_link_libraries(test_checkpoint PRIVATE
    ikafssn_sam_writer ikafssn_io ikafssn_core
    OpenSSL::Crypto hts
)
target_include_directories(test_checkpoint PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_SOURCE_DIR}/src"
)
add_test(NAME test_checkpoint COMMAND test_checkpoint)

# ikafssninfo integration test (requires NCBI Toolkit for index building)
add_executable(test_ikafssninfo test_ikafssninfo.cpp)
target_link_libraries(test_ikafssninfo PRIVATE
    ikafssn_builder ikafssn_index ikafssn_blastdb ikafssn_io ikafssn_util ikafssn_core
)
target_include_directories(test_ikafssninfo PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${NCBI_TOOLKIT_INCLUDE}"
    "${CMAKE_SOURCE_DIR}/src"
)
target_compile_definitions(test_ikafssninfo PRIVATE
    SOURCE_DIR="${CMAKE_SOURCE_DIR}"
)
add_test(NAME test_ikafssninfo COMMAND test_ikafssninfo)
