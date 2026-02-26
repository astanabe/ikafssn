function(add_ikafssn_test TEST_NAME TEST_SOURCE)
    add_executable(${TEST_NAME} ${TEST_SOURCE})
    target_link_libraries(${TEST_NAME} PRIVATE ikafssn_index ikafssn_io ikafssn_core)
    target_include_directories(${TEST_NAME} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endfunction()

add_ikafssn_test(test_varint test_varint.cpp)
add_ikafssn_test(test_kmer_encoding test_kmer_encoding.cpp)
add_ikafssn_test(test_ksx_io test_ksx_io.cpp)
add_ikafssn_test(test_kix_io test_kix_io.cpp)
add_ikafssn_test(test_kpx_io test_kpx_io.cpp)
add_ikafssn_test(test_khx_io test_khx_io.cpp)
target_link_libraries(test_khx_io PRIVATE ikafssn_util)

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
)
add_test(NAME test_builder COMMAND test_builder)

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
)
add_test(NAME test_packed_kmer_scanner COMMAND test_packed_kmer_scanner)

# --- Phase 3 tests ---

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
)
add_test(NAME test_stage1 COMMAND test_stage1)

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

# --- Phase 4 tests ---

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

# --- Phase 5 tests ---

# Protocol serialize/deserialize test (no external dependencies)
add_executable(test_protocol test_protocol.cpp)
target_link_libraries(test_protocol PRIVATE
    ikafssn_protocol ikafssn_util ikafssn_core
)
target_include_directories(test_protocol PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
add_test(NAME test_protocol COMMAND test_protocol)

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

# --- Phase 6 tests ---

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

# --- Phase 6 tests (retrieve) ---

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

# --- Phase 7 tests ---

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
