cmake_minimum_required(VERSION 3.31.6)
# set(CMAKE_EXPERIMENTAL_CXX_IMPORT_STD "0e5b6991-d74f-4b3d-a41c-cf096e0b2508")
# set(CMAKE_CXX_MODULE_STD cxx_std_23)
project(qchem VERSION 0.1.0 LANGUAGES CXX)


enable_language (Fortran)
enable_testing()
include(CTest)
include(GoogleTest)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/submodules/blaze/cmake )
include(Blaze_Import)

add_library(qchem)  
target_sources(qchem PRIVATE )
set_property(TARGET qchem PROPERTY CXX_STANDARD 20)
# set(CMAKE_CXX_STANDARD_REQUIRED YES)
# set(CMAKE_CXX_EXTENSIONS ON)

set(CXX_SCAN_FOR_MODULES on)

# Common compile flags for all targets
set(CXX_COMMON_COMPILE_OPTIONS
    $<$<CONFIG:Debug>:-g -O0>
    $<$<CONFIG:Release>:-O2 -pg>
)
set(CXX_COMMON_COMPILE_DEFINITIONS
    $<$<CONFIG:Debug>:-DDEBUG>
    $<$<CONFIG:Release>:-DNDEBUG>
)

target_compile_options(qchem PRIVATE ${CXX_COMMON_COMPILE_OPTIONS})
target_compile_definitions(qchem PUBLIC ${CXX_COMMON_COMPILE_DEFINITIONS})
target_link_libraries(qchem qcWaveFunction qcSCFAccelerator qcOrbitals qcHamiltonian qcChargeDensity qcFitting qcFactory qcMolecule_BS qcAtom_BS qcBasisSet qcLASolver qcCluster qcMesh qcSymmetry wignerSymbols lapack blas gtest )
add_subdirectory(src)
add_subdirectory(./UnitTests)
add_subdirectory(./gtkapp)

add_custom_target(allTests
    DEPENDS UTMain UTCommon UTMesh UTCluster
    COMMENT "Build all unit test executables"
)
# Disable embedded submodule tests so ctest only runs this project's tests.
set(BUILD_TESTING OFF CACHE BOOL "Disable embedded submodule tests" FORCE)
add_subdirectory(./submodules/wignerSymbols)
add_subdirectory(./submodules/libxc)
set(BUILD_TESTING ON CACHE BOOL "Enable testing for this project" FORCE)

Blaze_Import(
   QUIET
   BLAS on
   LAPACK on
   #THREADING Boost
   CACHE_SIZE auto
   VECTORIZATION on
   STORAGE_ORDER rowMajor
   THRESHOLD_DMATDVECMULT 100000UL
   THRESHOLD_SMP_DVECDVECADD 1000000UL
)

