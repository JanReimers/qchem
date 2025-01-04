cmake_minimum_required(VERSION 3.5.0)
project(qchem VERSION 0.1.0 LANGUAGES CXX)

file (GLOB_RECURSE qchem_SOURCES CONFIGURE_DEPENDS "src/*.C")
file (GLOB_RECURSE qchem_HEADERS CONFIGURE_DEPENDS "include/*.H")

# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/debug)
# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/release)
# set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/debug)
# set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/release)

add_library(qchem ${qchem_SOURCES}) 

target_compile_options(qchem PRIVATE -Wall -pg -fopenmp
    $<$<CONFIG:Debug>:-g -O0>
    $<$<CONFIG:Release>: -O2>
)
target_compile_definitions(qchem PRIVATE 
    $<$<CONFIG:Debug>:-DDEBUG>
    $<$<CONFIG:Release>:-DNDEBUG>
)
target_include_directories(qchem PRIVATE include ../oml/include)
target_link_libraries(qchem oml lapack gtest pthread wignerSymbols)
add_subdirectory(./UnitTests)

