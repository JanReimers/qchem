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
#      A panel case that FAILS today (an open tracker row) is listed in QCHEM_PANEL_RED below, with its row: it runs every time and PASSES only while gtest reports it failed -- the WILL_FAIL semantics (CLAUDE.md), but
#      as a PASS_REGULAR_EXPRESSION on gtest's own "[  FAILED  ] 1 test", so a red case whose filter matched nothing, or that
#      crashed, is red too (bare WILL_FAIL would invert "no tests were run" into a pass).  When one is fixed it goes red: remove
#      it from the list and close its row.
set(QCHEM_PANEL_RED
    "GPW_Na2Box.Γ_Pol_Null_Linear_SpinSeed_Converges_Panel"                # OpenWork §3 "Na₂ polarized singlet won't converge from AFM seed"
    "GPW_MnOBox.Γ_M6_GDM_Converges_Panel_Long"                # OpenWork §3 "GDM crawls at a near-degenerate frontier"
    "GPW_MnO2.Γ_U_Shub_Pol_Smear_Kerker_Pulay_Converges_Panel_Long" # OpenWork §3 "Kerker+Pulay overshoots β-MnO₂'s U = 1 eV basin"
    "GPW_MnO2.Γ_U_Shub_Pol_Smear_Kerker_ConvergesAtU2eV_Panel_Long"     # same row: control U = 2 eV limit-cycles under plain Kerker
    "ResponsePolarizability.GPW_Si_Γ_Pol_OneSpinAlpha_KeepsZeroMoment_Panel" # OpenWork §3 "Polarized Γ Si JUMPS" (panel #10)
)
set(QCHEM_PERIODIC
    "GPW_MnO.Γ_U_Shub_Pol_Smear_CP2K_Long"                                   # twin: GPW_MnO.Γ_U_Shub_Pol_Smear_Saved_CP2K_Long (vendored states)
    "ResponseLittleGroup.GPW_Si_k222_q222_Becke_ReducedEqualsFull_Long"     # twin: ResponseLittleGroup.GPW_Si_k222_q211_Becke_ReducedEqualsFull
    "ResponseLittleGroup.GPW_Si_k222_q211_ReducedEqualsFull_Long"           # (uniform raster) same twin: the cost is per q per solve
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
        TEST_FILTER -*_Long:*_Panel
        DISCOVERY_EXTRA_ARGS --gtest_output=
        PROPERTIES FAIL_REGULAR_EXPRESSION "no tests were run")
    gtest_discover_tests(${target}
        TEST_FILTER *_Long-*_Panel_Long:${_qchemPeriodic}
        DISCOVERY_EXTRA_ARGS --gtest_output=
        PROPERTIES FAIL_REGULAR_EXPRESSION "no tests were run" LABELS long TIMEOUT 3600)
    gtest_discover_tests(${target}
        TEST_FILTER ${_qchemPeriodic}
        DISCOVERY_EXTRA_ARGS --gtest_output=
        PROPERTIES FAIL_REGULAR_EXPRESSION "no tests were run" LABELS "long\;periodic" TIMEOUT 3600)
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
        PROPERTIES PASS_REGULAR_EXPRESSION "\\[  FAILED  \\] 1 test" LABELS "panel\;panelRed")
    gtest_discover_tests(${target}
        TEST_FILTER ${_qchemPanelRedLong}
        DISCOVERY_EXTRA_ARGS --gtest_output=
        PROPERTIES PASS_REGULAR_EXPRESSION "\\[  FAILED  \\] 1 test" LABELS "long\;panel\;panelRed" TIMEOUT 14400)
endfunction()
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

# QCHEM_ENZYME: the automatic-differentiation path (src/AD, qchem.AD; doc/GeometryOptimizationPlan.md, the AD spike).  Enzyme is
# TOOLCHAIN and OPTIONAL (controller ruling 2026-11-08): OFF builds everything, and every AD gradient REFUSES with NoGradient (an
# Outcome naming why); no hand-coded fallback derivative is kept beside the AD one.  ON needs an Enzyme built against THIS clang
# (LLVM 21): -DQCHEM_ENZYME=ON -DQCHEM_ENZYME_PLUGIN=<build>/Enzyme/ClangEnzyme-21.so -DQCHEM_ENZYME_INCLUDE=<src>/enzyme/include.
# The plugin rides ONLY qcAD's compile (it differentiates the kernels defined in its own TUs), so no other target's .pcm or
# compile command changes.
option(QCHEM_ENZYME "Build the Enzyme automatic-differentiation path (qchem.AD); OFF = AD gradients refuse with NoGradient" OFF)
set(QCHEM_ENZYME_PLUGIN  "" CACHE FILEPATH "Enzyme's clang plugin, ClangEnzyme-<LLVM major>.so (QCHEM_ENZYME=ON)")
set(QCHEM_ENZYME_INCLUDE "" CACHE PATH     "Enzyme's C++ headers, <Enzyme source>/enzyme/include (QCHEM_ENZYME=ON)")
if(QCHEM_ENZYME)
    if(NOT EXISTS "${QCHEM_ENZYME_PLUGIN}" OR NOT EXISTS "${QCHEM_ENZYME_INCLUDE}/enzyme/enzyme")
        message(FATAL_ERROR "QCHEM_ENZYME=ON needs QCHEM_ENZYME_PLUGIN (ClangEnzyme-21.so) and QCHEM_ENZYME_INCLUDE (holding enzyme/enzyme); "
                            "got '${QCHEM_ENZYME_PLUGIN}' and '${QCHEM_ENZYME_INCLUDE}'")
    endif()
    message(STATUS "QCHEM_ENZYME: AD on, plugin ${QCHEM_ENZYME_PLUGIN}")
endif()
# qchem_enzyme_sources(<target> <source>...): the TUs that REGISTER derivatives of their own leaf functions (doc/GeometryOptimizationPlan.md,
# the AD verdict: each Imp TU differentiates its own leaves in place).  Under QCHEM_ENZYME they compile with the plugin, its headers and
# the QCHEM_ENZYME macro, at -O2 -mprefer-vector-width=128 (after the tree's -O3, so they win): Enzyme 1a87ccc segfaulted at
# -O2 -march=native (256-bit vectors) and gave a wrong adjoint in a shape its -O3 produced (doc/Records/EnzymeSpike.md).  -march stays (a
# .pcm and its importers agree on target features); the flags cost the energy nothing measurable (the record §4 item 10).  Per FILE,
# so no other TU of the target changes.  OFF: nothing, and the files' #ifdef QCHEM_ENZYME branches refuse with NoGradient.
function(qchem_enzyme_sources target)
    if(QCHEM_ENZYME)
        set_source_files_properties(${ARGN} TARGET_DIRECTORY ${target} PROPERTIES
            COMPILE_OPTIONS "-fplugin=${QCHEM_ENZYME_PLUGIN};-O2;-mprefer-vector-width=128"
            COMPILE_DEFINITIONS QCHEM_ENZYME
            INCLUDE_DIRECTORIES "${QCHEM_ENZYME_INCLUDE}")
    endif()
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
if(QCHEM_RELCHECKED)
    set(CMAKE_CXX_FLAGS_RELEASE "${QCHEM_RELEASE_CODEGEN} -g -fno-omit-frame-pointer")
else()
    set(CMAKE_CXX_FLAGS_RELEASE "${QCHEM_RELEASE_CODEGEN} -g -fno-omit-frame-pointer -DNDEBUG")
endif()
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
target_link_libraries(qchem qcResponse qcWaveFunction qcSCFAccelerator qcOrbitals qcHamiltonian qcChargeDensity qcFitting qcFactory qcGaussian_BS qcRadial_BS qcLattice_BS qcPlaneWave_BS qcBasisSet qcLASolver qcStructure qcMesh qcSymmetry lapack blas)
add_subdirectory(src)
add_subdirectory(./IntegrationTests)
add_subdirectory(./CLIapps)

add_custom_target(allTests
    DEPENDS ITMain UTAD UTCommon UTCalculation UTHamiltonian UTChargeDensity UTPseudopotential UTMesh UTStructure UTSymmetry UTElConfig UTRadial_BS UTGaussian_BS UTLattice_BS UTLASolver UTSCFAccelerator UTSCFIterator UTResponse UTWaveFunction
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
# qcGaussian_BS is a shared library, so the static cint must be position-independent.
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
