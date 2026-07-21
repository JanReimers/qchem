cmake_minimum_required(VERSION 3.31.6)
# set(CMAKE_EXPERIMENTAL_CXX_IMPORT_STD "0e5b6991-d74f-4b3d-a41c-cf096e0b2508")
# set(CMAKE_CXX_MODULE_STD cxx_std_23)
project(qchem VERSION 0.1.0 LANGUAGES CXX)

if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    message(FATAL_ERROR
        "This project requires the Clang C++ toolchain. "
        "Use clang++/clang and flang, for example: "
        "CXX=/opt/LLVM-21.1.6-Linux-X64/bin/clang++ "
        "CC=/opt/LLVM-21.1.6-Linux-X64/bin/clang "
        "FC=/opt/LLVM-21.1.6-Linux-X64/bin/flang cmake -S . -B build/Release -G Ninja"
    )
endif()

set(BasisDataPath "BASISSET_DATA_PATH=\"${PROJECT_SOURCE_DIR}/BasisSetData\"" CACHE INTERNAL "BasisSet data path compile definition")
set(CommonDataPath "COMMON_DATA_PATH=\"${PROJECT_SOURCE_DIR}/src/Common/Data\"" CACHE INTERNAL "Common data path compile definition")
set(LatticeDataPath "LATTICE_DATA_PATH=\"${PROJECT_SOURCE_DIR}/src/BasisSet/Lattice_3D/Data\"" CACHE INTERNAL "Lattice_3D (GTH pseudopotential) data path compile definition")
set(PseudoDataPath "PSEUDO_DATA_PATH=\"${PROJECT_SOURCE_DIR}/src/Pseudopotential/Data\"" CACHE INTERNAL "Pseudopotential (GTH database) data path compile definition")
set(ChargeDensityDataPath "CHARGEDENSITY_DATA_PATH=\"${PROJECT_SOURCE_DIR}/src/ChargeDensity/Data\"" CACHE INTERNAL "ChargeDensity (SAD atomic-density database) data path compile definition")

enable_language (Fortran)
enable_testing()
include(CTest)
include(GoogleTest)
set(BUILD_GMOCK OFF CACHE BOOL "Disable GoogleMock for this project" FORCE)
add_subdirectory(./submodules/googletest)
link_directories(${PROJECT_BINARY_DIR}/lib)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/submodules/blaze/cmake )
include(Blaze_Import)

add_library(qchem)  
target_sources(qchem PRIVATE )
set_property(TARGET qchem PROPERTY CXX_STANDARD 20)
# set(CMAKE_CXX_STANDARD_REQUIRED YES)
# set(CMAKE_CXX_EXTENSIONS ON)

set(CXX_SCAN_FOR_MODULES on)

# RelChecked: a Release build (optimized -O2 codegen, i.e. the production object layout) but with
# asserts LIVE and blaze bounds checking on.  It catches layout-sensitive UB that -O0 Debug reshuffles
# away and that NDEBUG'd Release/ASan hide.  Configure with -DCMAKE_BUILD_TYPE=Release -DQCHEM_RELCHECKED=ON
# (Release, not a custom config name, because libxc's safeguard rejects unknown build types).
# See project_msym_layout_ub.
option(QCHEM_RELCHECKED "Release codegen with asserts live + blaze bounds checks (no NDEBUG)" OFF)

# Common compile flags for all targets
set(CXX_COMMON_COMPILE_OPTIONS
    $<$<CONFIG:Debug>:-g -O0>
    # Release keeps -O2; -g + frame pointers make it profilable (perf/callgrind) with negligible runtime
    # cost.  Dropped -pg (gprof): inaccurate on optimized code and superseded by perf/callgrind/PGO.
    $<$<CONFIG:Release>:-O2 -g -fno-omit-frame-pointer>
)
set(CXX_COMMON_COMPILE_DEFINITIONS
    $<$<CONFIG:Debug>:-DDEBUG>
)
if(QCHEM_RELCHECKED)
    # No -DNDEBUG (codebase asserts live) + blaze BLAZE_USER_ASSERT -> assert() (bounds checks live).
    list(APPEND CXX_COMMON_COMPILE_DEFINITIONS -DBLAZE_USER_ASSERTION=1)
else()
    list(APPEND CXX_COMMON_COMPILE_DEFINITIONS $<$<CONFIG:Release>:-DNDEBUG>)
endif()

target_compile_options(qchem PRIVATE ${CXX_COMMON_COMPILE_OPTIONS})
target_compile_definitions(qchem PUBLIC ${CXX_COMMON_COMPILE_DEFINITIONS})
target_link_libraries(qchem qcWaveFunction qcSCFAccelerator qcOrbitals qcHamiltonian qcChargeDensity qcFitting qcFactory qcMolecule_BS qcAtom_BS qcLattice_BS qcBasisSet qcLASolver qcStructure qcMesh qcSymmetry lapack blas)
add_subdirectory(src)
add_subdirectory(./UnitTests)

add_custom_target(allTests
    DEPENDS UTMain UTCommon UTChargeDensity UTPseudopotential UTMesh UTStructure UTSymmetry UTElConfig UTAtom_BS UTMolecule_BS UTLattice_3D_BS UTLASolver
    COMMENT "Build all unit test executables"
)
# Disable embedded submodule tests so ctest only runs this project's tests.
set(BUILD_TESTING OFF CACHE BOOL "Disable embedded submodule tests" FORCE)
add_subdirectory(./submodules/libxc)
# libcint (Apache-2.0): block-oriented Gaussian integral library, used by the PG_LibCint
# matrix-delivery evaluator.  Build it static with no Fortran interface (we call the C API directly).
set(ENABLE_STATIC ON CACHE BOOL "Build libcint as a static library" FORCE)
set(WITH_FORTRAN  OFF CACHE BOOL "libcint Fortran interface (unused)" FORCE)
add_subdirectory(./submodules/libcint)
# qcMolecule_BS is a shared library, so the static cint must be position-independent.
set_property(TARGET cint PROPERTY POSITION_INDEPENDENT_CODE ON)
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

# Python bindings (pybind/qchem_py) -- OFF by default so the normal build is untouched.
# Enable with: -DQCHEM_PYBIND=ON -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DPython_EXECUTABLE=<venv python>
option(QCHEM_PYBIND "Build the nanobind Python bindings (pybind/)" OFF)
if(QCHEM_PYBIND)
    add_subdirectory(pybind)
endif()

