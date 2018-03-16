cmake_minimum_required(VERSION 3.5)
project(read-filter)

set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH}
    "${CMAKE_CURRENT_SOURCE_DIR}/vendor/seqan/util/cmake")
set(SEQAN_INCLUDE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/vendor/seqan/include/")

# Dependencies
find_package(ZLIB)
find_package(BZip2)
find_package(OpenMP)
find_package(SeqAn REQUIRED)

# Setup SeqAn
include_directories(${SEQAN_INCLUDE_DIRS})
add_definitions(${SEQAN_DEFINITIONS})
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

add_executable(read-filter
    src/main.cpp
)

target_link_libraries(read-filter ${SEQAN_LIBRARIES})
