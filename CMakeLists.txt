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

add_subdirectory(src)

# -Wno-unknown-pragmas
target_compile_options(qchem PRIVATE 
    $<$<CONFIG:Debug>:-g -O0>
    $<$<CONFIG:Release>: -O2 >
)
target_compile_definitions(qchem PRIVATE 
    $<$<CONFIG:Debug>:-DDEBUG>
    $<$<CONFIG:Release>:-DNDEBUG>
)
target_link_libraries(qchem qcWaveFunction qcSCFAccelerator qcOrbitals qcHamiltonian qcChargeDensity qcFitting qcFactory qcMolecule_BS qcAtom_BS qcBasisSet qcLASolver_blaze qcCluster qcMesh qcSymmetry omlm wignerSymbols lapack blas gtest )
add_subdirectory(./UnitTests)
add_subdirectory(./gtkapp)
add_subdirectory(./submodules/oml)
add_subdirectory(./submodules/wignerSymbols)

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

