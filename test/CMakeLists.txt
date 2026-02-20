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
