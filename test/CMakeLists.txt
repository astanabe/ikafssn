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
