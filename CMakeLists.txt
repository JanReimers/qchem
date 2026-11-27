cmake_minimum_required(VERSION 3.31.6)
# `import std` (user ruling 2026-11-25, with libc++ below): every target is C++23 (a C++20 TU cannot load a C++23 BMI) and may import
# std / std.compat as modules.  CMake still gates this behind an EXPERIMENTAL UUID that CHANGES with each CMake release: the one below is
# CMake 4.2's.  On another CMake configure stops here; to upgrade, take the new UUID from that CMake's Help/dev/experimental.rst.
if(NOT CMAKE_VERSION VERSION_GREATER_EQUAL 4.2 OR CMAKE_VERSION VERSION_GREATER_EQUAL 4.3)
    message(FATAL_ERROR "import std is gated by a per-release CMake UUID; this tree holds CMake 4.2's and this is CMake ${CMAKE_VERSION}.\n"
        "-- To fix this, please configure with CMake 4.2.x, or update CMAKE_EXPERIMENTAL_CXX_IMPORT_STD in the root CMakeLists.txt to this "
        "CMake's UUID (Help/dev/experimental.rst) and widen this check.")
endif()
set(CMAKE_EXPERIMENTAL_CXX_IMPORT_STD "d0edc3af-4c50-42ea-a356-e2862fe7a444")
set(CMAKE_CXX_MODULE_STD ON)
set(CMAKE_CXX_STANDARD 23)
# ONE INPUT: the C++ compiler (user, 2026-11-23: "external users who try and build, just need to provide a compiler path ... and our
# cmake files handle the rest").  The C and Fortran compilers default to clang and flang BESIDE the stated clang++ (the dependencies
# build C, and Fortran is enabled below), so a tree never mixes the stated LLVM with the system cc.  A stated CC / FC / CMAKE_*_COMPILER wins.
set(_qchemCxx "${CMAKE_CXX_COMPILER}")
if(NOT _qchemCxx)
    set(_qchemCxx "$ENV{CXX}")
endif()
if(_qchemCxx MATCHES "clang\\+\\+(-[0-9]+)?$")
    get_filename_component(_qchemLLVMBin "${_qchemCxx}" DIRECTORY)
    if(NOT DEFINED CMAKE_C_COMPILER AND NOT DEFINED ENV{CC} AND EXISTS "${_qchemLLVMBin}/clang")
        set(CMAKE_C_COMPILER "${_qchemLLVMBin}/clang" CACHE FILEPATH "C compiler: clang beside the stated clang++")
    endif()
    if(NOT DEFINED CMAKE_Fortran_COMPILER AND NOT DEFINED ENV{FC} AND EXISTS "${_qchemLLVMBin}/flang")
        set(CMAKE_Fortran_COMPILER "${_qchemLLVMBin}/flang" CACHE FILEPATH "Fortran compiler: flang beside the stated clang++")
    endif()
endif()
# THE C++ STANDARD LIBRARY IS LLVM's libc++ (user ruling 2026-11-25, "Go ahead with libc++"; lane G study laneG-import-std-feasibility.md).
# It is part of the toolchain, so it is the project's choice and never an input: it goes into the INITIAL CMAKE_CXX_FLAGS (every C++
# compile AND link, try_compile and the dependencies included).  It is chosen before project() so the compiler checks run with it.  The
# check after project() refuses a tree whose cached flags predate this (a reconfigured libstdc++ tree): libc++ and libstdc++ objects
# must never mix in one binary.
string(APPEND CMAKE_CXX_FLAGS_INIT " -stdlib=libc++")
project(qchem VERSION 0.1.0 LANGUAGES CXX)
# Release unless stated: the tests, the anchors and the timings are all Release numbers (RelChecked is Release + QCHEM_RELCHECKED).
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type (default Release)" FORCE)
    message(STATUS "CMAKE_BUILD_TYPE not stated: Release")
endif()

if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    message(FATAL_ERROR
        "This project requires the Clang C++ toolchain. "
        "State clang++ alone (C and Fortran are derived from it), on a fresh tree, for example: "
        "cmake -S . -B build/Release -G Ninja -DCMAKE_CXX_COMPILER=/opt/LLVM-23.1.4-Linux-X64/bin/clang++"
    )
endif()
# FULL BMIs.  Clang 23 writes REDUCED BMIs by default, and those do not merge the duplicate Blaze partial specializations that
# qchem.Types and qchem.Blaze both include textually in their global-module fragments: MultTraitEval2, SolveTraitEval2 ... turn
# ambiguous and 59 objects fail (lane G, 2026-11-23).  -fno-modules-reduced-bmi restores the full BMIs clang 21 writes by default;
# clang 21 does not know the flag, so it is added only where the compiler accepts it.  A compile option (not CMAKE_CXX_FLAGS),
# so it stays off the link line.
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-fno-modules-reduced-bmi QCHEM_CXX_HAS_NO_REDUCED_BMI)
if(QCHEM_CXX_HAS_NO_REDUCED_BMI)
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-fno-modules-reduced-bmi>)
endif()
# libc++: the tree must really be on it, and every binary must LOAD the libc++ of this clang (the /opt tarball keeps it in its per-target
# runtime dir, which is not on the loader's default search, so its directory joins CMAKE_BUILD_RPATH here, before any target and before
# the OpenMP runtime check below).
include(CheckCXXSourceCompiles)
check_cxx_source_compiles("#include <version>\n#ifndef _LIBCPP_VERSION\n#error not libc++\n#endif\nint main(){}" QCHEM_ON_LIBCXX)
if(NOT QCHEM_ON_LIBCXX)
    message(FATAL_ERROR "This project builds on LLVM's libc++ (user ruling 2026-11-25), and this tree compiles against another standard "
        "library (CMAKE_CXX_FLAGS='${CMAKE_CXX_FLAGS}' has no -stdlib=libc++: the tree was configured before the switch, or the flags were "
        "stated by hand).\n-- To fix this, please configure a FRESH (wiped) tree, or add -stdlib=libc++ to a stated CMAKE_CXX_FLAGS.")
endif()
execute_process(COMMAND ${CMAKE_CXX_COMPILER} -stdlib=libc++ -print-file-name=libc++.so.1 OUTPUT_VARIABLE _qchemLibcxx OUTPUT_STRIP_TRAILING_WHITESPACE)
if(NOT IS_ABSOLUTE "${_qchemLibcxx}" OR NOT EXISTS "${_qchemLibcxx}")
    message(FATAL_ERROR "libc++: ${CMAKE_CXX_COMPILER} does not ship libc++.so.1 (-print-file-name answered '${_qchemLibcxx}').\n"
        "-- To fix this, please state a clang++ whose LLVM install includes libc++ (the /opt/LLVM-23.1.4-Linux-X64 tarball does).")
endif()
get_filename_component(_qchemLibcxx "${_qchemLibcxx}" REALPATH)
get_filename_component(_qchemLibcxxDir "${_qchemLibcxx}" DIRECTORY)
list(APPEND CMAKE_BUILD_RPATH "${_qchemLibcxxDir}")
message(STATUS "C++ standard library: libc++ (${_qchemLibcxxDir}, in the build RPATH)")
# The Python bindings are a shared module, so every object linked into it is position independent: QCHEM_PYBIND implies PIC
# (without it the thread_local qchem::detail::g_outerDepth takes a local-exec TLS relocation no shared object can use).
# Declared here, before the first target, because CMAKE_POSITION_INDEPENDENT_CODE is read when a target is created.
option(QCHEM_PYBIND "Build the nanobind Python bindings (pybind/)" OFF)
if(QCHEM_PYBIND)
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

# OpenMP: used ONLY for the GPW collocate/integrate pair-loop parallelism in qcGaussian_BS
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
    # The /opt LLVM 23.1.4 tarball (the project's toolchain since 2026-11-24) ships its own omp.h and libomp, and find_package finds
    # them unaided.  HISTORY: the /opt LLVM 21.1.6 tarball lacked omp.h, and its trees pointed find_package at apt's libomp-dev:
    #   cmake -DOpenMP_CXX_FLAGS="-fopenmp -I/usr/lib/llvm-21/lib/clang/21/include" \
    #         -DOpenMP_CXX_LIB_NAMES=omp -DOpenMP_omp_LIBRARY=/usr/lib/x86_64-linux-gnu/libomp.so .
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    add_compile_definitions(QCHEM_OPENMP)
    link_libraries(${OpenMP_omp_LIBRARY})
    # The libomp a binary LOADS must be the one it LINKS (2026-11-23, lane G, clang 23: the tree linked /opt's libomp by full path, CMake
    # left that directory out of RUNPATH -- it is one of the /opt clang's implicit link dirs -- and every binary silently ran on apt's older
    # /usr/lib libomp).  cmake/OpenMPRuntimeCheck.cpp is linked the tree's way and prints the libomp it loads at run time.  When the loader's
    # default search already finds the linked one (apt's libomp in /usr/lib) nothing is added: a system dir at the head of RUNPATH would
    # reorder every other library's lookup.  Otherwise the library's directory joins CMAKE_BUILD_RPATH, and the check must then pass.
    separate_arguments(_qchemOmpFlags UNIX_COMMAND "${OpenMP_CXX_FLAGS}")
    get_filename_component(_qchemOmpLinked "${OpenMP_omp_LIBRARY}" REALPATH)
    get_filename_component(_qchemOmpDir "${OpenMP_omp_LIBRARY}" DIRECTORY)
    function(qchem_omp_loaded out rpaths)
        set(_rpathOpts "")
        foreach(_d IN LISTS rpaths)
            list(APPEND _rpathOpts "-Wl,-rpath,${_d}")
        endforeach()
        try_run(_run _compile
                SOURCES "${PROJECT_SOURCE_DIR}/cmake/OpenMPRuntimeCheck.cpp" NO_CACHE
                COMPILE_DEFINITIONS ${_qchemOmpFlags}
                LINK_OPTIONS ${_qchemOmpFlags} ${_rpathOpts}
                LINK_LIBRARIES "${OpenMP_omp_LIBRARY}"
                COMPILE_OUTPUT_VARIABLE _compileOut RUN_OUTPUT_STDOUT_VARIABLE _loaded)  # stdout only: a TSan libomp WARNS on stderr
        if(NOT _compile OR NOT _run EQUAL 0)
            message(FATAL_ERROR "OpenMP: a program linked against ${OpenMP_omp_LIBRARY} with '${OpenMP_CXX_FLAGS}' does not build or run:\n"
                "${_compileOut}\n-- To fix this, please state -DOpenMP_omp_LIBRARY=<the libomp.so of this clang> (and OpenMP_CXX_FLAGS) on a fresh tree.")
        endif()
        set(${out} "${_loaded}" PARENT_SCOPE)
    endfunction()
    qchem_omp_loaded(_qchemOmpLoaded "${CMAKE_BUILD_RPATH}")
    set(_qchemOmpHow "by the loader's default search")
    if(NOT _qchemOmpLoaded STREQUAL _qchemOmpLinked)
        list(APPEND CMAKE_BUILD_RPATH "${_qchemOmpDir}")
        qchem_omp_loaded(_qchemOmpLoaded "${CMAKE_BUILD_RPATH}")
        set(_qchemOmpHow "through the build RPATH ${_qchemOmpDir}")
    endif()
    if(NOT _qchemOmpLoaded STREQUAL _qchemOmpLinked)
        message(FATAL_ERROR "OpenMP: the tree links ${_qchemOmpLinked}, but a program built the same way LOADS ${_qchemOmpLoaded} at run time, "
            "even with ${_qchemOmpDir} in its RPATH (LD_LIBRARY_PATH='$ENV{LD_LIBRARY_PATH}', which overrides RPATH).  Every binary would run "
            "silently on another OpenMP runtime.\n-- To fix this, please unset LD_LIBRARY_PATH for the build, or state "
            "-DOpenMP_omp_LIBRARY=<the libomp.so you mean to run on> on a fresh tree.")
    endif()
    message(STATUS "OpenMP runtime: ${_qchemOmpLoaded}, loaded as linked (${_qchemOmpHow})")
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

# Blaze -> BLAS for the dense products.  ⚠ DEFAULT **ON** -- see the "DEFAULT ON since 2026-08-15" note
# below.  (This line used to read "OFF by default; -DQCHEM_BLAZE_BLAS=ON to try it", which stayed behind
# when the option flipped and then misled doc/Records/ParallelAndOraclePlan.md 1.2 into planning a session around
# turning on something that had been on for three weeks.  Corrected 2026-09-06.)
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
set(LatticeDataPath "LATTICE_DATA_PATH=\"${PROJECT_SOURCE_DIR}/src/BasisSet/Lattice/Data\"" CACHE INTERNAL "Lattice_3D (GTH pseudopotential) data path compile definition")
set(PseudoDataPath "PSEUDO_DATA_PATH=\"${PROJECT_SOURCE_DIR}/src/Pseudopotential/Data\"" CACHE INTERNAL "Pseudopotential (GTH database) data path compile definition")
set(StructureDataPath "STRUCTURE_DATA_PATH=\"${PROJECT_SOURCE_DIR}/src/Structure/Data\"" CACHE INTERNAL "Pre-defined structures (molecules/materials JSON) data path compile definition")
set(ChargeDensityDataPath "CHARGEDENSITY_DATA_PATH=\"${PROJECT_SOURCE_DIR}/src/ChargeDensity/Data\"" CACHE INTERNAL "ChargeDensity (SAD atomic-density database) data path compile definition")

enable_language (Fortran)
enable_testing()
# The written module rule as a test: no PUBLIC module unit may re-export an .Internal. module (V1.20c/d).
add_test(NAME InternalReexportAudit COMMAND ${PROJECT_SOURCE_DIR}/scripts/audit-internal-reexports ${PROJECT_SOURCE_DIR}/src)
# The tier-1 abstract hubs are PROTECTED faces (V-INTERFACE-TIERS step 3): each one's exported declarations and re-exports against
# the checked-in snapshot; a change fails until the snapshot is regenerated (--update) in a commit citing the design ruling.
add_test(NAME HubInterfaceAudit COMMAND ${PROJECT_SOURCE_DIR}/scripts/audit-hub-interfaces ${PROJECT_SOURCE_DIR}/src
                                        ${PROJECT_SOURCE_DIR}/scripts/hub-interfaces.snapshot)
# Every environment variable the library reads must be on the settings reference page (src/Common/Settings.dox) -- D-ENV.
add_test(NAME SettingsDocAudit COMMAND ${PROJECT_SOURCE_DIR}/scripts/audit-settings-doc ${PROJECT_SOURCE_DIR}/src)
add_test(NAME PPChoiceAudit COMMAND ${PROJECT_SOURCE_DIR}/scripts/audit-pp-choice ${PROJECT_SOURCE_DIR}/src)
# D-ENV 6e: no environment variable outside the resource/diagnostics allowlist may be read anywhere (the deck is the only way to set a number)
add_test(NAME GetenvAudit COMMAND ${PROJECT_SOURCE_DIR}/scripts/audit-getenv ${PROJECT_SOURCE_DIR})
# The BasisSet G-tag rule as a test (doc/Records/BasisSetTaxonomyPlan.md §1.6 rule 2, step 3): no qchem.BasisSet.Gaussian.Point.*
# module imports a lattice; no qchem.BasisSet.Radial.* module imports a lattice or the molecular point-group symmetry.
add_test(NAME BasisSetGTagAudit COMMAND ${PROJECT_SOURCE_DIR}/scripts/audit-basisset-gtags ${PROJECT_SOURCE_DIR}/src/BasisSet)
include(CTest)
include(GoogleTest)
# qchem_discover_tests(<exe>): gtest_discover_tests with the two guards every test exe in this tree needs
# (2026-09-21, after ctest silently skipped 48 Γ-named SCF tests for five days):
#   1. DISCOVERY_EXTRA_ARGS --gtest_output=  -- TEXT-listing discovery.  CMake >= 3.30 lists through
#      --gtest_output=json and googletest escapes every non-ASCII BYTE of a name as \u00XX, so a non-ASCII
#      name comes back double-encoded and its filter matches nothing.  The trailing empty flag suppresses the
#      JSON file and the module parses the plain listing, bytes intact.
#   2. FAIL_REGULAR_EXPRESSION "no tests were run" -- googletest exits 0 when a filter matches NOTHING (it
#      only WARNS: "filter ... did not match any test; no tests were run"), which defeats the point of a test
#      runner; ctest turns that line into a FAILURE, so the next encoding or renaming accident is red, not green.
# The real fix -- a runner that cannot pass an empty selection -- is doc/ModuleToolchainPlan.md's gtest item.
#   3. THE `_Long` TIER (D-LONG, ruled 2026-10-02): a test whose name ENDS in `_Long` is over the 60 s CPU budget.
#      It is discovered by a second pass and labelled `long`, so the regression sweep is `ctest -j8 -LE long` and
#      the long gates run as `ctest -j8 -L long`.  (A plain `ctest` still runs both tiers; ctest has no
#      per-test default-off.)  The two filters are disjoint, so the ctest count is the same as one pass.
#   4. THE CONVERGENCE PANEL (doc/SCFRobustnessPlan.md §4, ruled 2026-10-29): a name ending `_Panel` is a fast panel case (label
#      `panel`, still in the regular sweep); `_Panel_Long` is the full tier (labels `long;panel`, so `-LE long` skips it).
#      `ctest -L panel` runs the panel -- ON DEMAND, under the heavy lock, and not in the routine long tier (user ruling
#      2026-10-30): the routine long tier is `ctest -j3 -L long -LE "panel|periodic"` (item 5), and a landing adds `-L panel` only when it
#      touches the SCF path (SCFAccelerator / SCFIterator / mixers / Calculation's SCF).  The `_Panel_Long` passes carry TIMEOUT
#      14400 (4 h): serial `GPW_MnO2.Γ_U_Shub_Pol_Smear_Kerker_Converges_Panel_Long` outran 7200 s (4268 s at 4 threads).  The
#      plain `_Long` pass carries TIMEOUT 3600: its slowest case, `GPW_MnO.Γ_U_Shub_Pol_Smear_CP2K_Long`, took 1294-1347 s serial,
#      too close to ctest's 1500 s default -- an explicit margin until the cheap twins replace these (OpenWork §3 "Every `_Long`
#      test re-proves its claim at the cheapest resolution").
#   5. THE PERIODIC TIER (OpenWork §3 "Every `_Long` test re-proves its claim at the cheapest resolution", controller ruling
#      2026-10-30): a COLD `_Long` gate whose claim a cheap twin now carries -- an energy anchor resumed from a vendored state
#      (decks/states/), an identity at a smaller q-mesh -- is listed in QCHEM_PERIODIC below beside its twin.  It keeps the
#      claim honest end to end and is the only thing that REFRESHES a vendored state (by an explicit cp, never as a sweep side
#      effect), so it is labelled `long;periodic` and runs with the panel: `ctest -L "panel|periodic"`.  The routine long tier
#      is `ctest -j3 -L long -LE "panel|periodic"`.
#      Its TIMEOUT is 14400 (the panel's 4 h): the cold NaF relaxation takes 8659 s alone on 4 threads, rutile 3320 s (lane B 2026-11-21).
#      A panel case that FAILS today (an open tracker row) is listed in QCHEM_PANEL_RED below, with its row: it runs every time and PASSES only while gtest reports it failed -- the WILL_FAIL semantics (CLAUDE.md), but
#      as a PASS_REGULAR_EXPRESSION on gtest's own "[  FAILED  ] 1 test", so a red case whose filter matched nothing, or that
#      crashed, is red too (bare WILL_FAIL would invert "no tests were run" into a pass).  When one is fixed it goes red: remove
#      it from the list and close its row.
#      ⚠ A RED CASE PASSES ON AN ASSERTION FAILURE, NEVER ON A THROW (lane B 2026-11-09): both red tiers also carry
#      FAIL_REGULAR_EXPRESSION "C++ exception with description" (ctest's FAIL regex wins over PASS), so a case that throws for an
#      unrelated reason -- a deck refusal, a rotted setup -- reads Failed, not Passed.  Earned by GPW_MnO's ordering gate: disabled
#      since 2026-10-04, it had rotted (R2.5 made it state siteSpins) and threw in 43 ms, which the PASS regex alone read as Passed.
set(QCHEM_PANEL_RED
    "GPW_Na2Box.Γ_Pol_Null_Linear_SpinSeed_Converges_Panel"                # OpenWork §3 "Na₂ polarized singlet won't converge from AFM seed"
    "GPW_MnOBox.Γ_M6_GDM_Converges_Panel_Long"                # OpenWork §3 "GDM crawls at a near-degenerate frontier"
    "GPW_MnO2.Γ_U_Shub_Pol_Smear_Kerker_Pulay_Converges_Panel_Long" # OpenWork §3 "Kerker+Pulay overshoots β-MnO₂'s U = 1 eV basin"
    "GPW_MnO2.Γ_U_Shub_Pol_Smear_Kerker_ConvergesAtU2eV_Panel_Long"     # same row: control U = 2 eV limit-cycles under plain Kerker
    "ResponsePolarizability.GPW_Si_Γ_Pol_OneSpinAlpha_KeepsZeroMoment_Panel" # OpenWork §3 "Polarized Γ Si JUMPS" (panel #10)
)
#   6. KNOWN RED (user ruling 2026-11-08, CLAUDE.md DISABLED_/WILL_FAIL rule): a NON-panel case that fails today for a reason an
#      open OpenWork §3 row names, and whose gate must NOT be loosened, is listed in QCHEM_KNOWN_RED with its row.  It runs in its
#      own tier as the panel-red cases do (PASSES only on gtest's "[  FAILED  ] 1 test", never on a throw), label `knownRed` (a `_Long` one also `long`): `ctest -L knownRed`
#      lists them.  A sweep report lists them separately with their reason.  When one goes red (fixed), remove it and close the row.
set(QCHEM_KNOWN_RED
    "GPW_TiO2.k222_CP2KsCellIsStationary_Long"                 # OpenWork §3 "TiO2 k222 energy vs CP2K is MESH-limited on the stated (40,29) Becke mesh": the ENERGY line fails (-142 uHa), stationarity passes
)
list(JOIN QCHEM_KNOWN_RED ":" _qchemKnownRed)
set(_qchemKnownRedFast ${QCHEM_KNOWN_RED})        # a _Long known-red case keeps the long label and timeout
list(FILTER _qchemKnownRedFast EXCLUDE REGEX "_Long$")
set(_qchemKnownRedLong ${QCHEM_KNOWN_RED})
list(FILTER _qchemKnownRedLong INCLUDE REGEX "_Long$")
list(JOIN _qchemKnownRedFast ":" _qchemKnownRedFast)
list(JOIN _qchemKnownRedLong ":" _qchemKnownRedLong)
set(QCHEM_PERIODIC
    "GPW_MnO.Γ_U_Shub_Pol_Smear_CP2K_Long"                                   # twin: GPW_MnO.Γ_U_Shub_Pol_Smear_Saved_CP2K_Long (vendored states)
    "GPW_MnO.k222_Shub_Pol_Smear_CP2K_Long"                                  # cold k-mesh oracle (~2 h serial); twin: GPW_Mn2Box.k211_Becke_Pol_Smear_EqualsItsSupercellAtGamma_Long (the k-phase identity)
    "ResponseLittleGroup.GPW_Si_k222_q222_Becke_ReducedEqualsFull_Long"     # twin: ResponseLittleGroup.GPW_Si_k222_q211_Becke_ReducedEqualsFull
    "ResponseLittleGroup.GPW_Si_k222_q211_ReducedEqualsFull_Long"           # (uniform raster) same twin: the cost is per q per solve
    "GPW_TiO2.k222_RelaxesToCP2KsCell_Long"                                 # cold rutile relaxation (55 min, 4 threads); twin: GPW_TiO2.k222_CP2KsCellIsStationary_Long (one point)
    "GPW_NaF.k222_RelaxesToCP2KsCell_Long"                                  # cold NaF relaxation (2.4 h, 4 threads); twin: GPW_NaF.k222_CP2KsCellIsStationary_Long (one point)
    "GPW_NaF.k222_CP2KsCellIsStationary_Long"                               # NOT a cheap twin: ONE point is 2328 s alone (OpenWork §4, the Becke force), past 3600 s beside two others
    "GPW_Si.k333_RelaxesToCP2KsCell_Long"                                   # cold Si relaxation on the DEFAULT accuracy (33 min, the pilot); twin: GPW_Si.k333_CP2KsCellIsStationary_Long
    "GPW_AlN.k222_CP2KsCellIsStationary_Long"                               # NOT a cheap twin: ONE point on stated grids is 3205 s on 4 threads (lane B 2026-11-22), 89% of the long TIMEOUT; no cold k222 gate (OpenWork §3 AlN row)
)
list(JOIN QCHEM_PERIODIC ":" _qchemPeriodic)
list(JOIN QCHEM_PANEL_RED ":" _qchemPanelRed)
set(_qchemPanelRedFast ${QCHEM_PANEL_RED})
list(FILTER _qchemPanelRedFast EXCLUDE REGEX "_Panel_Long$")
set(_qchemPanelRedLong ${QCHEM_PANEL_RED})
list(FILTER _qchemPanelRedLong INCLUDE REGEX "_Panel_Long$")
list(JOIN _qchemPanelRedFast ":" _qchemPanelRedFast)
list(JOIN _qchemPanelRedLong ":" _qchemPanelRedLong)
function(qchem_discover_tests target)
    gtest_discover_tests(${target}
        TEST_FILTER -*_Long:*_Panel:${_qchemKnownRed}
        DISCOVERY_EXTRA_ARGS --gtest_output=
        PROPERTIES FAIL_REGULAR_EXPRESSION "no tests were run")
    gtest_discover_tests(${target}
        TEST_FILTER *_Long-*_Panel_Long:${_qchemPeriodic}:${_qchemKnownRed}
        DISCOVERY_EXTRA_ARGS --gtest_output=
        PROPERTIES FAIL_REGULAR_EXPRESSION "no tests were run" LABELS long TIMEOUT 3600)
    gtest_discover_tests(${target}
        TEST_FILTER ${_qchemPeriodic}
        DISCOVERY_EXTRA_ARGS --gtest_output=
        PROPERTIES FAIL_REGULAR_EXPRESSION "no tests were run" LABELS "long\;periodic" TIMEOUT 14400)
    gtest_discover_tests(${target}
        TEST_FILTER *_Panel-${_qchemPanelRed}
        DISCOVERY_EXTRA_ARGS --gtest_output=
        PROPERTIES FAIL_REGULAR_EXPRESSION "no tests were run" LABELS panel)
    gtest_discover_tests(${target}
        TEST_FILTER *_Panel_Long-${_qchemPanelRed}
        DISCOVERY_EXTRA_ARGS --gtest_output=
        PROPERTIES FAIL_REGULAR_EXPRESSION "no tests were run" LABELS "long\;panel" TIMEOUT 14400)
    gtest_discover_tests(${target}
        TEST_FILTER ${_qchemPanelRedFast}
        DISCOVERY_EXTRA_ARGS --gtest_output=
        PROPERTIES PASS_REGULAR_EXPRESSION "\\[  FAILED  \\] 1 test" FAIL_REGULAR_EXPRESSION "C\\+\\+ exception with description" LABELS "panel\;panelRed")
    gtest_discover_tests(${target}
        TEST_FILTER ${_qchemPanelRedLong}
        DISCOVERY_EXTRA_ARGS --gtest_output=
        PROPERTIES PASS_REGULAR_EXPRESSION "\\[  FAILED  \\] 1 test" FAIL_REGULAR_EXPRESSION "C\\+\\+ exception with description" LABELS "long\;panel\;panelRed" TIMEOUT 14400)
    if(NOT "${_qchemKnownRedFast}" STREQUAL "")   # an EMPTY filter would discover every test as known-red
        gtest_discover_tests(${target}
            TEST_FILTER ${_qchemKnownRedFast}
            DISCOVERY_EXTRA_ARGS --gtest_output=
            PROPERTIES PASS_REGULAR_EXPRESSION "\\[  FAILED  \\] 1 test" FAIL_REGULAR_EXPRESSION "C\\+\\+ exception with description" LABELS knownRed)
    endif()
    if(NOT "${_qchemKnownRedLong}" STREQUAL "")
        gtest_discover_tests(${target}
            TEST_FILTER ${_qchemKnownRedLong}
            DISCOVERY_EXTRA_ARGS --gtest_output=
            PROPERTIES PASS_REGULAR_EXPRESSION "\\[  FAILED  \\] 1 test" FAIL_REGULAR_EXPRESSION "C\\+\\+ exception with description" LABELS "long\;knownRed" TIMEOUT 3600)
    endif()
endfunction()
# THE SOURCE DEPENDENCIES: one pin table (cmake/Dependencies.cmake, doc/SharedDependenciesPlan.md); submodules/ is no longer read (step 3 deletes it).
include(${PROJECT_SOURCE_DIR}/cmake/Dependencies.cmake)
link_directories(${PROJECT_BINARY_DIR}/lib)
list(APPEND CMAKE_MODULE_PATH ${QCHEM_BLAZE_DIR}/cmake )
include(Blaze_Import)

add_library(qchem)  
target_sources(qchem PRIVATE )
set_property(TARGET qchem PROPERTY CXX_STANDARD 23)
# set(CMAKE_CXX_STANDARD_REQUIRED YES)
# set(CMAKE_CXX_EXTENSIONS ON)

set(CXX_SCAN_FOR_MODULES on)

# RelChecked: a Release build (optimized -O2 codegen, i.e. the production object layout) but with
# asserts LIVE and blaze bounds checking on.  It catches layout-sensitive UB that -O0 Debug reshuffles
# away and that NDEBUG'd Release/ASan hide.  Configure with -DCMAKE_BUILD_TYPE=Release -DQCHEM_RELCHECKED=ON
# (Release, not a custom config name, because libxc's safeguard rejects unknown build types).
# See project_msym_layout_ub.
option(QCHEM_RELCHECKED "Release codegen with asserts live + blaze bounds checks (no NDEBUG)" OFF)

# ENZYME IS REQUIRED (user ruling 2026-11-09: "Yes I think we make it required. Seems like a small burden for a lot of
# functionality.").  Enzyme (an LLVM plugin) differentiates the IR of the TUs that ask for a derivative (src/AD, and the per-file
# registrations below; doc/GeometryOptimizationPlan.md, the AD spike; doc/Records/EnzymeSpike.md).  There is no build without it:
# configure STOPS with a reason when the plugin is missing, built against another LLVM, older than the pin, or does not differentiate.
#   QCHEM_ENZYME_PLUGIN   ClangEnzyme-<LLVM major>.so, built against THIS clang (found under ~/Code/Enzyme-<pin>/build/Enzyme or
#                         ~/Code/Enzyme-<pin>/build-llvm<LLVM major>/Enzyme by default: one Enzyme checkout, one build dir per LLVM)
#   QCHEM_ENZYME_INCLUDE  <Enzyme source>/enzyme/include (found under ~/Code/Enzyme-<pin> by default)
# The default checkout is NAMED FOR THE PIN (controller ruling 2026-11-13): a pin move builds Enzyme in a new ~/Code/Enzyme-<pin> beside the
# old one, so no tree ever compiles against a plugin swapped under it; a tree whose cache still names the old checkout fails the pin
# check below LOUDLY, and `-UQCHEM_ENZYME_PLUGIN -UQCHEM_ENZYME_INCLUDE` on its reconfigure finds the new one.
# PIN: Enzyme cdd0d1c2 (v0.0.309-23, 2026-11-12; moved from 871de19f on 2026-11-13).  It holds the fixes of both upstream defects
# the user filed: EnzymeAD/Enzyme#3747 (a silently wrong reverse-mode gradient: a callee filling a loop-local array through an inner
# loop with a data-dependent trip count; fixed by cdd0d1c2 itself, #3810) and #3748 (a segfault under 256-bit vectors; fixed by
# d0cc6f5f, #3754).  Both filed repros are exact on it at -O1..-O3, with and without -march=native (doc/Records/EnzymeSpike.md §6).
# Configure checks the include tree's git checkout holds the pin.  Every Enzyme derivative keeps its finite-difference gate all the same.
# The 1a87ccc-era per-file workaround flags (-O2 -mprefer-vector-width=128) are gone with the defects (controller ruling 2026-11-13).
if(DEFINED CACHE{QCHEM_ENZYME})
    unset(QCHEM_ENZYME CACHE)
    message(STATUS "QCHEM_ENZYME is retired and ignored: Enzyme is required (user ruling 2026-11-09)")
endif()
set(QCHEM_ENZYME_PIN cdd0d1c2)
if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    message(FATAL_ERROR "Enzyme is required and is a clang plugin; this compiler is ${CMAKE_CXX_COMPILER_ID}")
endif()
string(REGEX MATCH "^[0-9]+" _qchemLLVMMajor "${CMAKE_CXX_COMPILER_VERSION}")
foreach(_v QCHEM_ENZYME_PLUGIN QCHEM_ENZYME_INCLUDE)   # an empty cached value (a tree configured before the ruling) searches again
    if(DEFINED CACHE{${_v}} AND "${${_v}}" STREQUAL "")
        unset(${_v} CACHE)
    endif()
endforeach()
find_file(QCHEM_ENZYME_PLUGIN ClangEnzyme-${_qchemLLVMMajor}.so
          HINTS "$ENV{HOME}/Code/Enzyme-${QCHEM_ENZYME_PIN}/build/Enzyme" "$ENV{HOME}/Code/Enzyme-${QCHEM_ENZYME_PIN}/build-llvm${_qchemLLVMMajor}/Enzyme"
          DOC "Enzyme's clang plugin, ClangEnzyme-<LLVM major>.so, built against this clang (required)")
find_path(QCHEM_ENZYME_INCLUDE enzyme/enzyme HINTS "$ENV{HOME}/Code/Enzyme-${QCHEM_ENZYME_PIN}/enzyme/include"
          DOC "Enzyme's C++ headers, <Enzyme source>/enzyme/include (required)")
set(_qchemEnzymeHow "Build Enzyme ${QCHEM_ENZYME_PIN} against this clang (LLVM ${_qchemLLVMMajor}) in $ENV{HOME}/Code/Enzyme-${QCHEM_ENZYME_PIN} "
    "(the build recipe: doc/Records/EnzymeSpike.md) and reconfigure with\n"
    "  -UQCHEM_ENZYME_PLUGIN -UQCHEM_ENZYME_INCLUDE   (the default finds it there), or name another checkout with\n"
    "  -DQCHEM_ENZYME_PLUGIN=<Enzyme>/build/Enzyme/ClangEnzyme-${_qchemLLVMMajor}.so -DQCHEM_ENZYME_INCLUDE=<Enzyme>/enzyme/include\n")
if(NOT EXISTS "${QCHEM_ENZYME_PLUGIN}" OR NOT EXISTS "${QCHEM_ENZYME_INCLUDE}/enzyme/enzyme")
    message(FATAL_ERROR "Enzyme is REQUIRED (user ruling 2026-11-09) and was not found:\n"
        "  QCHEM_ENZYME_PLUGIN  = '${QCHEM_ENZYME_PLUGIN}' (want ClangEnzyme-${_qchemLLVMMajor}.so)\n"
        "  QCHEM_ENZYME_INCLUDE = '${QCHEM_ENZYME_INCLUDE}' (want the directory holding enzyme/enzyme)\n" ${_qchemEnzymeHow})
endif()
# The pin: the include tree's checkout must hold QCHEM_ENZYME_PIN (Enzyme's own version number is not maintained in its CMake).
get_filename_component(_qchemEnzymeSrc "${QCHEM_ENZYME_INCLUDE}/../.." ABSOLUTE)
find_package(Git QUIET)
if(GIT_FOUND)
    execute_process(COMMAND "${GIT_EXECUTABLE}" -C "${_qchemEnzymeSrc}" merge-base --is-ancestor ${QCHEM_ENZYME_PIN} HEAD
                    RESULT_VARIABLE _qchemPinRC OUTPUT_QUIET ERROR_QUIET)
    execute_process(COMMAND "${GIT_EXECUTABLE}" -C "${_qchemEnzymeSrc}" describe --tags --always
                    OUTPUT_VARIABLE _qchemEnzymeRev OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    if(NOT _qchemPinRC EQUAL 0)
        message(FATAL_ERROR "Enzyme at ${_qchemEnzymeSrc} ('${_qchemEnzymeRev}') does not hold the pinned commit ${QCHEM_ENZYME_PIN} "
            "(older Enzymes silently give a wrong reverse-mode gradient, EnzymeAD/Enzyme#3747, or segfault under AVX2, #3748).\n" ${_qchemEnzymeHow})
    endif()
else()
    set(_qchemEnzymeRev "unverified: no git")
    message(WARNING "Cannot check the Enzyme pin ${QCHEM_ENZYME_PIN} (no git); the smoke test below still runs")
endif()
# The smoke test: the plugin must load into this clang and differentiate (a plugin built against another LLVM fails here, not mid-build).
# Its -O2 is the smoke program's own (Enzyme runs in the optimisation pipeline), not a flag of the tree's TUs.
# The smoke program RUNS, so it carries the tree's build RPATH (libc++, libomp): try_run does not apply CMAKE_BUILD_RPATH itself.
set(_qchemEnzymeRpath "")
foreach(_d IN LISTS CMAKE_BUILD_RPATH)
    list(APPEND _qchemEnzymeRpath "-Wl,-rpath,${_d}")
endforeach()
try_run(_qchemEnzymeRun _qchemEnzymeCompile
        SOURCES "${PROJECT_SOURCE_DIR}/cmake/EnzymeSmoke.cpp"
        COMPILE_DEFINITIONS "-fplugin=${QCHEM_ENZYME_PLUGIN}" -O2
        LINK_OPTIONS ${_qchemEnzymeRpath}
        CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${QCHEM_ENZYME_INCLUDE}"
        COMPILE_OUTPUT_VARIABLE _qchemEnzymeCompileOut RUN_OUTPUT_VARIABLE _qchemEnzymeRunOut)
if(NOT _qchemEnzymeCompile OR NOT _qchemEnzymeRun EQUAL 0)
    message(FATAL_ERROR "The Enzyme plugin ${QCHEM_ENZYME_PLUGIN} does not differentiate under this clang (LLVM ${_qchemLLVMMajor}): "
        "d(x^2)/dx at 3 gave '${_qchemEnzymeRunOut}' (want 6).\n${_qchemEnzymeCompileOut}\n" ${_qchemEnzymeHow})
endif()
message(STATUS "Enzyme: ${_qchemEnzymeRev} (pin ${QCHEM_ENZYME_PIN}), plugin ${QCHEM_ENZYME_PLUGIN}; smoke d(x^2)/dx = 6")
# qchem_enzyme_sources(<target> <source>...): the TUs that REGISTER derivatives of their own leaf functions (doc/GeometryOptimizationPlan.md,
# the AD verdict: each Imp TU differentiates its own leaves in place).  They compile with the plugin and its headers, at the TREE's codegen
# flags: the per-file -O2 -mprefer-vector-width=128 they carried until 2026-11-13 dodged EnzymeAD/Enzyme#3747 (a wrong adjoint at -O3) and
# #3748 (an AVX2 segfault), both fixed at the pin.  Per FILE, so no other TU of the target changes.
# OBJECT_DEPENDS on the plugin: ninja does not see a -fplugin= .so as an input, so without it a rebuilt plugin at the same path left these
# objects as the OLD plugin compiled them (controller ruling 2026-11-13).  src/AD registers all of qcAD's TUs through it.
function(qchem_enzyme_sources target)
    set_source_files_properties(${ARGN} TARGET_DIRECTORY ${target} PROPERTIES
        COMPILE_OPTIONS "-fplugin=${QCHEM_ENZYME_PLUGIN}"
        INCLUDE_DIRECTORIES "${QCHEM_ENZYME_INCLUDE}"
        OBJECT_DEPENDS "${QCHEM_ENZYME_PLUGIN}")
endfunction()

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
# -O from its peers by accident and was invisible in perf call graphs); qcRadial_BS kept a private -O2 copy
# missing -g.  Adding -march=native then FAILED TO BUILD outright -- "precompiled file qchem.Math.pcm was
# compiled with the target feature '+avx2' but the current translation unit is not" -- because a .pcm and
# every TU consuming it must agree on target features.  In CMAKE_CXX_FLAGS_<CONFIG> a target CANNOT opt
# out, and changing them rewrites every compile command, so the whole cosmos rebuilds.
#
# QCHEM_ARCH_EXPERIMENT picks the Release codegen.
#
# ★ NATIVE IS THE DEFAULT AS OF 2026-09-04.  The earlier reading (193.3 / 191.1 / 181.2 s single-shot,
# with an O3 repeat at 197.3) could not separate ~6% from ~3-4% whole-run noise, so it was left off.
# Re-measured on all three MnO rows with doc/Benchmark.md's OWN PRINTED COMMAND (VA span, annealed,
# serial) -- CPU seconds, O3 -> NATIVE:
#
#     ALL DEFAULTS      427.1 -> 386.1 s   -9.6%
#     QCHEM_BECKE_XC=0  169.9 -> 155.7 s   -8.4%
#     CP2K_COMPAT probe 357.9 -> 338.9 s   -5.3%
#
# and E_tot is BIT-IDENTICAL on every row (-61.40297551 / -61.40358773, both matching the 08-28 bank).
# -O3 vs -O2 remains INSIDE the noise and is still not a measured win.
#
# ⚠ AN EARLIER 09-04 READING OF "-5.2%" WAS TAKEN ON THE WRONG CONFIGURATION (default basis, single SCF
# stage -- see doc/Benchmark.md's retraction).  Same direction, wrong system; these numbers replace it.
#
# ⚠ ON THE TWO OBJECTIONS THAT KEPT IT OFF:
#   * "results move" (-march licenses FMA CONTRACTION: a*b+c rounds ONCE, not twice).  Real, but Etot came
#     out BIT-IDENTICAL across all six runs above (-61.41280709).  What it did move: three ROUNDOFF FLOORS
#     in the BSpline/DBCache tests -- checks on quantities that are mathematically EXACTLY zero, whose
#     thresholds had been fitted 1.05-1.5x above one build's rounding.  Those are widened to a few ulps
#     with the reason recorded at each site.  No physics anchor was re-pinned, and none needed to be.
#   * "the Phi bucket is ~15% worse" -- DEAD.  ALL DEFAULTS *is* the Becke/Phi route, and it is the row
#     NATIVE wins by the MOST (-9.6%).  Whatever that 2026-08-20 regression was, the profile has been
#     rewritten twice since and it no longer shows at the run level.
#
# ⚠ SCOPE: only C++ gets it (CMAKE_CXX_FLAGS_<CONFIG>).  Verified on the configured tree: 527 C++ TUs carry
# `-O3 -march=native`, while 374 .c and 1 .f TU keep plain -O3 (their languages are untouched by design)
# and googletest's 2 .cc files opt out in their own subproject -- harmless, since none of them import a
# .pcm and so none can hit the "compiled with target feature '+avx2'" mismatch this file warns about.
# (Binary PORTABILITY is not a concern here -- R&D project, no binary consumers.)
# ⚠ AN EXISTING BUILD TREE KEEPS ITS OLD VALUE.  `set(... CACHE ...)` only supplies a default when the
# entry does not exist yet, so a tree configured before 2026-09-04 stays on whatever it cached (O3) and
# this default changes NOTHING for it.  To move one over:
#     cmake -DQCHEM_ARCH_EXPERIMENT=NATIVE build/Release      # or delete the entry / the tree
# Check what a tree ACTUALLY compiles with -- never what you think it should:
#     grep -m1 QCHEM_ARCH_EXPERIMENT: build/Release/CMakeCache.txt
# and confirm the flags reach every TU with exactly ONE -O (the 2026-08-20 drift this file records):
#     python3 -c "import json,re,collections;db=json.load(open('build/Release/compile_commands.json'));\
#     print(collections.Counter(tuple(re.findall(r'(?<= )-O\w*',e.get('command') or ' '.join(e['arguments']))) for e in db))"
# A clean tree prints {('-O3',): N} -- one -O per TU, no duplicates.
set(QCHEM_ARCH_EXPERIMENT "NATIVE" CACHE STRING "Release codegen: NATIVE (default; -O3 -march=native, licenses FMA contraction) | O3 | O2")
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
# -DNDEBUG here too (it also reaches the submodules' C++), but NOT under QCHEM_RELCHECKED: an unconditional one silently
# killed every assert in the "asserts live" tree (found 2026-10-31 setting up the daily RelChecked sweep).
# RelChecked also HARDENS (user ruling 2026-11-26: "I think it makes sense to turn on _LIBCPP_HARDENING_MODE_EXTENSIVE and
# -fsanitize=array-bounds for the 2:30am debug runs"):
#  - libc++'s EXTENSIVE hardening checks std::array / std::vector / std::span / string_view element access (operator[], front, back);
#  - -fsanitize=array-bounds checks a fixed-size C array INDEXED BY A RUNTIME VALUE, including one INSIDE a class: an overrun that stays
#    inside the object, which ASan cannot see (the 2026-11-03 Hermite3 LMAX overrun, which corrupted Li q3's local PP by 9-30 mHa);
#  - -fsanitize-trap=array-bounds makes a failed check a TRAP (SIGILL at the faulting line, readable in a core or gdb): no runtime library.
# Both are in the GLOBAL flags, like the codegen flags above: a macro and a sanitizer are part of a BMI's configuration, so every TU of a
# tree must agree.  Release stays unhardened (its numbers and timings are the anchors).
set(QCHEM_RELCHECKED_HARDENING "-D_LIBCPP_HARDENING_MODE=_LIBCPP_HARDENING_MODE_EXTENSIVE -fsanitize=array-bounds -fsanitize-trap=array-bounds")
if(QCHEM_RELCHECKED)
    set(CMAKE_CXX_FLAGS_RELEASE "${QCHEM_RELEASE_CODEGEN} -g -fno-omit-frame-pointer ${QCHEM_RELCHECKED_HARDENING}")
else()
    set(CMAKE_CXX_FLAGS_RELEASE "${QCHEM_RELEASE_CODEGEN} -g -fno-omit-frame-pointer -DNDEBUG")
endif()
set(CMAKE_CXX_FLAGS_DEBUG   "-g -O0")

# A class that marks SOME overrides `override` and leaves others unmarked is an ERROR (user ruling 2026-11-19, lane H's
# warning cleanup).  clang raises -Winconsistent-missing-override only then, and in a module interface it is re-reported by
# every importing TU (956 of the 957 warnings in the b414d40a landing build, from 146 lines).  An override that exists on ONE
# template instantiation only lives in a conditional CRTP mixin (Composite_HFSystem, HamiltonianImp_RealBlock,
# PeriodicIrrepCD_RealEnergy), never as an unmarked member that overrides nothing on the other.  C++ only.
add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-Werror=inconsistent-missing-override>)

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
target_link_libraries(qchem qcOptimize qcResponse qcWaveFunction qcSCFAccelerator qcOrbitals qcHamiltonian qcChargeDensity qcFitting qcFactory qcGaussian_BS qcRadial_BS qcLattice_BS qcPlaneWave_BS qcBasisSet qcLASolver qcStructure qcMesh qcSymmetry lapack blas)
add_subdirectory(src)
add_subdirectory(./IntegrationTests)
add_subdirectory(./CLIapps)

# The CLI apps ride every sweep build: a stale rundeck/valgen against a rebuilt libqchem segfaulted (lane G) or silently dropped
# its report (lane H, 2026-11-11), because no sweep or nightly relinked them.  A new app in CLIapps/ joins this list.
add_custom_target(allApps
    COMMENT "Build all CLI apps"
)
add_dependencies(allApps scfrun valgen gpwprobe rundeck gth2upf)   # add_dependencies, not DEPENDS: the executable audit below reads it
# Benchmarks and profiling drivers: not tests, but they must COMPILE on every landing (lane B 2026-11-13: MnDBench had been broken since
# c926486ca, 2026-11-04, because nothing in allTests built it -- the _NOT_BUILT trap's silent-rot class).
add_custom_target(allBenches
    COMMENT "Build all benchmark / profiling drivers"
)
add_dependencies(allBenches MnDBench)   # add_dependencies, not DEPENDS: the executable audit below reads it
add_custom_target(allTests
    COMMENT "Build all unit test executables"
)
add_dependencies(allTests allApps allBenches ITMain UTAD UTCommon UTOptimize UTCalculation UTHamiltonian UTChargeDensity UTPseudopotential UTMesh UTStructure UTSymmetry UTElConfig UTRadial_BS UTGaussian_BS UTLattice_BS UTLASolver UTSCFAccelerator UTSCFIterator UTResponse UTWaveFunction)   # add_dependencies, not DEPENDS: the executable audit below reads it
# libxc and libcint (cmake/Dependencies.cmake), after the project's own targets, their tests off.
qchem_build_c_dependencies()

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
# Enable with: -DQCHEM_PYBIND=ON [-DPython_EXECUTABLE=<venv python>, to bind for that Python]; the option and its PIC sit near the top.
if(QCHEM_PYBIND)
    add_subdirectory(pybind)
endif()


# THE LINKER ENFORCES LIBRARY DEPENDENCIES (OpenWork "-Wl,-z,defs", 2026-11-01): every project SHARED library links with -z defs, so
# a call into a library it does not DECLARE (target_link_libraries) fails at ITS OWN link -- by default a .so may leave symbols
# undefined and the omission surfaces only, if at all, when some executable links.  The targets are walked from the project's own
# trees (the root, src/, IntegrationTests/, CLIapps/), so a new library is covered with no list to keep; submodules and pybind/ are
# not touched.  A sanitizer build (its runtime resolves only in the executable) would have to drop it; this tree has none.
function(qchem_shared_libs_under dir out)
    get_property(targets DIRECTORY ${dir} PROPERTY BUILDSYSTEM_TARGETS)
    get_property(subdirs DIRECTORY ${dir} PROPERTY SUBDIRECTORIES)
    set(found "")
    foreach(t IN LISTS targets)
        get_target_property(type ${t} TYPE)
        if(type STREQUAL "SHARED_LIBRARY")
            list(APPEND found ${t})
        endif()
    endforeach()
    foreach(d IN LISTS subdirs)
        qchem_shared_libs_under(${d} sub)
        list(APPEND found ${sub})
    endforeach()
    set(${out} ${found} PARENT_SCOPE)
endfunction()
get_property(qchem_root_targets DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY BUILDSYSTEM_TARGETS)
set(qchem_defs_libs "")
foreach(t IN LISTS qchem_root_targets)
    get_target_property(type ${t} TYPE)
    if(type STREQUAL "SHARED_LIBRARY")
        list(APPEND qchem_defs_libs ${t})
    endif()
endforeach()
foreach(d src IntegrationTests CLIapps)
    qchem_shared_libs_under(${CMAKE_CURRENT_SOURCE_DIR}/${d} libs)
    list(APPEND qchem_defs_libs ${libs})
endforeach()
foreach(t IN LISTS qchem_defs_libs)
    target_link_options(${t} PRIVATE "LINKER:-z,defs")
endforeach()

# EVERY PROJECT EXECUTABLE IS BUILT BY allTests (lane B 2026-11-13).  An executable no sweep builds rots silently: MnDBench broke on
# 2026-11-04 and nobody saw it for nine days, the class the _NOT_BUILT trap belongs to (CLAUDE.md).  So configure walks the project's own
# trees (the root, src/, IntegrationTests/, CLIapps/; submodules and pybind/ are not touched) for EXECUTABLE targets and STOPS, naming
# them, when one is not reachable from allTests through the add_dependencies of the custom targets (allApps, allBenches, ...).  Fix: add it to the
# matching list above.
function(qchem_executables_under dir out)
    get_property(targets DIRECTORY ${dir} PROPERTY BUILDSYSTEM_TARGETS)
    get_property(subdirs DIRECTORY ${dir} PROPERTY SUBDIRECTORIES)
    set(found "")
    foreach(t IN LISTS targets)
        get_target_property(type ${t} TYPE)
        if(type STREQUAL "EXECUTABLE")
            list(APPEND found ${t})
        endif()
    endforeach()
    foreach(d IN LISTS subdirs)
        qchem_executables_under(${d} sub)
        list(APPEND found ${sub})
    endforeach()
    set(${out} ${found} PARENT_SCOPE)
endfunction()
set(qchem_all_exes "")
foreach(t IN LISTS qchem_root_targets)
    get_target_property(type ${t} TYPE)
    if(type STREQUAL "EXECUTABLE")
        list(APPEND qchem_all_exes ${t})
    endif()
endforeach()
foreach(d src IntegrationTests CLIapps)
    qchem_executables_under(${CMAKE_CURRENT_SOURCE_DIR}/${d} exes)
    list(APPEND qchem_all_exes ${exes})
endforeach()
set(qchem_reached "")
set(qchem_frontier allTests)
while(qchem_frontier)
    list(POP_FRONT qchem_frontier t)
    if(t IN_LIST qchem_reached OR NOT TARGET ${t})
        continue()
    endif()
    list(APPEND qchem_reached ${t})
    get_target_property(deps ${t} MANUALLY_ADDED_DEPENDENCIES)
    if(deps)
        list(APPEND qchem_frontier ${deps})
    endif()
endwhile()
set(qchem_unbuilt "")
foreach(e IN LISTS qchem_all_exes)
    if(NOT e IN_LIST qchem_reached)
        list(APPEND qchem_unbuilt ${e})
    endif()
endforeach()
if(qchem_unbuilt)
    message(FATAL_ERROR "These executables are not built by allTests, so no sweep or landing would ever compile them: ${qchem_unbuilt}.  "
        "Add each to the add_dependencies list of allApps, allBenches or allTests in the root CMakeLists.txt (the rule beside this check).")
endif()
