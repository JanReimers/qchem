cmake_minimum_required(VERSION 3.5.0)
project(qchem VERSION 0.1.0 LANGUAGES CXX)

enable_language (Fortran)
file (GLOB_RECURSE qchem_SOURCES CONFIGURE_DEPENDS "src/*.C")
file (GLOB_RECURSE qchem_FORTRAN CONFIGURE_DEPENDS "src/*.f")
file (GLOB_RECURSE qchem_HEADERS CONFIGURE_DEPENDS "include/*.H")
file (GLOB_RECURSE qchem_TEMPLATES CONFIGURE_DEPENDS "src/BasisSet/Atom/radial/BSpline/*.Ci")

# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/debug)
# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/release)
# set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/debug)
# set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/release)

add_library(qchem) 
target_sources(qchem PRIVATE ${qchem_SOURCES} ${qchem_FORTRAN})

# -Wno-unknown-pragmas
target_compile_options(qchem PRIVATE -Wall -Wno-unknown-pragmas
    $<$<CONFIG:Debug>:-g -O0>
    $<$<CONFIG:Release>: -O2>
)
target_compile_definitions(qchem PRIVATE 
    $<$<CONFIG:Debug>:-DDEBUG>
    $<$<CONFIG:Release>:-DNDEBUG>
)
target_include_directories(qchem PRIVATE include ../oml/include BSplinebasis/include src/BasisSet/Atom/radial/BSpline)
target_link_libraries(qchem oml lapack gtest wignerSymbols)
add_subdirectory(./UnitTests)
add_subdirectory(./gtkapp)

