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

# OpenMP: used ONLY for the GPW collocate/integrate pair-loop parallelism in qcMolecule_BS
# (PG_Cart_MnD::NR_Evaluator).  Opt-in at runtime (GPW_OMP_THREADS>1) and serial by default, so it never
# changes the bit-anchors.  The code gates on our own QCHEM_OPENMP macro (set alongside the flag), a
# holdover from the -fopenmp=libgomp era (link-compat only, no _OPENMP) that remains the project-owned
# gate now that real -fopenmp defines _OPENMP everywhere.
find_package(OpenMP QUIET)
if(OpenMP_CXX_FOUND)
    # OpenMP must be TRULY GLOBAL under C++20 modules: the OpenMP language mode is part of a PCM's
    # configuration and clang refuses to MIX (importing a non-OpenMP PCM into an -fopenmp TU errors
    # "OpenMP support differs in precompiled file", and vice versa) -- so per-target scoping cannot
    # build, and a per-target option VAR is not applied uniformly enough.  (The -O/-march codegen flags
    # below learned the SAME lesson the hard way in 2026-08-20 and now live in CMAKE_CXX_FLAGS_<CONFIG>
    # for exactly this reason -- target features are part of a PCM's configuration too.)
    # Hence: the flag goes into CMAKE_CXX_FLAGS (every C++ TU, Fortran/C untouched), the define into
    # add_compile_definitions (directory-inherited everywhere), and the runtime into link_libraries
    # (every target declared after this point).  QCHEM_OPENMP gates the (opt-in via GPW_OMP_THREADS,
    # serial-by-default) parallel code.
    # On a toolchain whose clang lacks omp.h (the /opt tarball), point find_package at apt's libomp-dev:
    #   cmake -DOpenMP_CXX_FLAGS="-fopenmp -I/usr/lib/llvm-21/lib/clang/21/include" \
    #         -DOpenMP_CXX_LIB_NAMES=omp -DOpenMP_omp_LIBRARY=/usr/lib/x86_64-linux-gnu/libomp.so .
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    add_compile_definitions(QCHEM_OPENMP)
    link_libraries(${OpenMP_omp_LIBRARY})
    # Blaze SMP must stay OFF: real -fopenmp defines _OPENMP in every TU, which flips Blaze's
    # shared-memory parallelization ON by default -- an unintended side effect of the global OpenMP
    # mode.  Blaze then wraps EVERY smpAssign (even tiny, below-threshold ones) in BLAZE_PARALLEL_SECTION,
    # whose guard flag is a GLOBAL non-atomic bool -- two of our GPW pair-loop workers doing any Blaze
    # vector copy race on it and one throws std::runtime_error("Nested parallel sections detected"),
    # which escapes the OpenMP region as std::terminate (the intermittent threaded-NaF abort,
    # doc/GPWPlan1.md).  Serial Blaze is also what every committed anchor was calibrated against, and
    # the parallelism strategy is ONE flat OpenMP region at the pair loop -- never nested regions.
    add_compile_definitions(BLAZE_USE_SHARED_MEMORY_PARALLELIZATION=0)
endif()

# Blaze -> BLAS for the dense products (OFF by default; -DQCHEM_BLAZE_BLAS=ON to try it).
#
# Blaze's own kernels are SIMD but build at the SSE2 baseline here (BUILD_MARCH_NATIVE is OFF), and the
# GPW hot path is now dominated by big dense products -- the XC-mesh Phi GEMMs, O(npts n^2) per SCF
# iteration (doc/GPWPlan1.md "THE RUNTIME GAP, MEASURED").  With OpenBLAS installed AND pinned to one
# thread (qchem::PinBlasToOneThread), routing those into zgemm/dgemm composes correctly with our own
# threading: OUR OpenMP region provides the parallelism, BLAS provides a fast serial kernel per block --
# one level of parallelism, deterministic reductions.
#
# It is a global define, not the Blaze cmake target's INTERFACE one: our TUs reach Blaze through
# `import qchem.Blaze` and never LINK the Blaze target, so the interface define does not reach them.
# BLAZE_BLAS_IS_PARALLEL stays 0 (its default) -- that is the pin, restated where Blaze can see it.
#
# DEFAULT ON since 2026-08-15: measured 34.1 vs 2.55 GFlop/s on the XC-mesh product, and `ctest -j8`
# is 716/716 green with it (540 s vs 615 s).  A BLAS kernel sums in a different order than Blaze's, so
# last-ULP values move -- that is the only cost, and the sweep says the anchors absorb it.  Turn it off
# with -DQCHEM_BLAZE_BLAS=OFF on a box without OpenBLAS or to reproduce pre-2026-08-15 digits exactly.
# NB the code CHOOSES ITS SHAPE from this: with BLAS mode on, PWTerms' quadrature does ONE whole-matrix
# product (Blaze only dispatches when the DESTINATION is a plain matrix -- see the note there); with it
# off, the hand-blocked threaded triangular form.  So flipping this is not purely a numerics switch.
option(QCHEM_BLAZE_BLAS "Route Blaze's dense products through the system BLAS (needs OpenBLAS)" ON)
if(QCHEM_BLAZE_BLAS)
    add_compile_definitions(BLAZE_BLAS_MODE=1 BLAZE_BLAS_IS_PARALLEL=0)
    link_libraries(openblas)      # every target declared after this point resolves cblas_* itself
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

# QCHEM_ARCH_EXPERIMENT: the -O3 / -march=native A/B (doc/OpenWork.md Step 3), OPT-IN and OFF by default.
# Release's -O2 is a DELIBERATE choice (below), not a flag-order accident, so the bar for changing it is a
# measurement -- and -march=native is not free: it licenses FMA CONTRACTION, so floating-point results
# MOVE (a*b+c fused rounds once, not twice).  Any bit-identical anchor must be re-checked, not assumed.
# It also pins the binary to this box's ISA.  Configure with -DQCHEM_ARCH_EXPERIMENT=ON to measure.
# ---- CODEGEN / ABI FLAGS: ONE PLACE, NO OPT-OUT (user, 2026-08-20) ------------------------------------
# These go into CMAKE_CXX_FLAGS_<CONFIG>, NOT a per-target option variable, for the reason the OpenMP flag
# above already records: "per-target option VARS are not applied uniformly enough".  That was not a style
# preference -- it is a CORRECTNESS requirement under C++20 modules, and the old per-target
# CXX_COMMON_COMPILE_OPTIONS had already drifted in three ways before anything forced the issue
# (measured 2026-08-20): qcMesh and Atom/Evaluators applied it NOWHERE and inherited a bare -O3 with no
# frame pointers (so the library owning BeckeCutoff, the hottest function in the code, was on a different
# -O from its peers by accident and was invisible in perf call graphs); qcAtom_BS kept a private -O2 copy
# missing -g.  Adding -march=native then FAILED TO BUILD outright -- "precompiled file qchem.Math.pcm was
# compiled with the target feature '+avx2' but the current translation unit is not" -- because a .pcm and
# every TU consuming it must agree on target features.  In CMAKE_CXX_FLAGS_<CONFIG> a target CANNOT opt
# out, and changing them rewrites every compile command, so the whole cosmos rebuilds.
#
# QCHEM_ARCH_EXPERIMENT picks the Release codegen.  MEASURED on the MnO fast loop (serial, doc/OpenWork.md
# Step 3): the old mixed state 193.3 s, O3 191.1 s, NATIVE 181.2 s -- BUT a repeat of the O3 build read
# 197.3 s, so single-sample whole-run timings carry ~3-4% noise on this box.  READ THAT AS: -O3 vs -O2 is
# INSIDE the noise and is NOT a measured win; -march's ~6% sits just outside it and would need repeats to
# stand up.  The default is O3 because unification forces ONE level on everyone and O3 is the only choice
# that DEMOTES no target (qcMesh was already there) -- not because it was measured faster.
# NATIVE is NOT the default: -march licenses FMA CONTRACTION, so results move, and it makes the Phi bucket
# ~15% WORSE (19.2 -> 22.1 s).  It belongs in the batched anchor re-pin alongside V1.22 / FittingCleanup K,
# not on its own.  (Binary PORTABILITY is not a concern here -- R&D project, no binary consumers.)
set(QCHEM_ARCH_EXPERIMENT "O3" CACHE STRING "Release codegen: O3 (default) | O2 | NATIVE (-O3 -march=native; moves FP results)")
set_property(CACHE QCHEM_ARCH_EXPERIMENT PROPERTY STRINGS O3 O2 NATIVE)
if(QCHEM_ARCH_EXPERIMENT STREQUAL "O2")
    set(QCHEM_RELEASE_CODEGEN "-O2")
elseif(QCHEM_ARCH_EXPERIMENT STREQUAL "NATIVE")
    set(QCHEM_RELEASE_CODEGEN "-O3 -march=native")
else()
    set(QCHEM_RELEASE_CODEGEN "-O3")
endif()
# -g + frame pointers keep Release PROFILABLE (perf/callgrind) at negligible cost -- that is what made the
# Step 3 cost attribution possible at all.  Dropped -pg (gprof): inaccurate on optimized code.
set(CMAKE_CXX_FLAGS_RELEASE "${QCHEM_RELEASE_CODEGEN} -g -fno-omit-frame-pointer -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG   "-g -O0")

set(CXX_COMMON_COMPILE_DEFINITIONS
    $<$<CONFIG:Debug>:-DDEBUG>
)
if(QCHEM_RELCHECKED)
    # No -DNDEBUG (codebase asserts live) + blaze BLAZE_USER_ASSERT -> assert() (bounds checks live).
    list(APPEND CXX_COMMON_COMPILE_DEFINITIONS -DBLAZE_USER_ASSERTION=1)
else()
    list(APPEND CXX_COMMON_COMPILE_DEFINITIONS $<$<CONFIG:Release>:-DNDEBUG>)
endif()

# (OpenMP handling lives next to find_package(OpenMP) near the top: C++20 modules force it into
#  CMAKE_CXX_FLAGS / add_compile_definitions / link_libraries -- truly global, not per-target vars.)

target_compile_definitions(qchem PUBLIC ${CXX_COMMON_COMPILE_DEFINITIONS})
target_link_libraries(qchem qcWaveFunction qcSCFAccelerator qcOrbitals qcHamiltonian qcChargeDensity qcFitting qcFactory qcMolecule_BS qcAtom_BS qcLattice_BS qcBasisSet qcLASolver qcStructure qcMesh qcSymmetry lapack blas)
add_subdirectory(src)
add_subdirectory(./IntegrationTests)
add_subdirectory(./CLIapps)

add_custom_target(allTests
    DEPENDS ITMain UTCommon UTHamiltonian UTChargeDensity UTPseudopotential UTMesh UTStructure UTSymmetry UTElConfig UTAtom_BS UTMolecule_BS UTLattice_3D_BS UTLASolver UTSCFAccelerator UTSCFIterator
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

