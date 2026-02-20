# Core library (header-only for now, but define as INTERFACE target)
add_library(ikafssn_core INTERFACE)
target_include_directories(ikafssn_core INTERFACE
    "${CMAKE_CURRENT_SOURCE_DIR}"
)

# Index library
add_library(ikafssn_index STATIC
    index/ksx_writer.cpp
    index/ksx_reader.cpp
    index/kix_writer.cpp
    index/kix_reader.cpp
    index/kpx_writer.cpp
    index/kpx_reader.cpp
)
target_link_libraries(ikafssn_index PUBLIC ikafssn_core ikafssn_io)

# IO library
add_library(ikafssn_io STATIC
    io/mmap_file.cpp
)
target_link_libraries(ikafssn_io PUBLIC ikafssn_core)
