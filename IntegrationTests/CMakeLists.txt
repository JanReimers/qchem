
link_directories(
    $<$<CONFIG:Debug>:  submodules/libxc/build>
    $<$<CONFIG:Release>:  submodules/libxc/build>
)

# THE GPW HARNESS (doc/TestSuitePlan.md phase 2): the recipe vocabulary + instruments the GPW integration
# tests and CLIapps/gpwprobe share.  Its own library so both compile one copy; no gtest in it.
add_library(qcGPW_Harness)
target_sources(qcGPW_Harness PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES GPW/Harness.C)
target_compile_definitions(qcGPW_Harness PUBLIC ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET qcGPW_Harness PROPERTY CXX_STANDARD 20)
target_include_directories(qcGPW_Harness PRIVATE
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    ${PROJECT_SOURCE_DIR}/submodules/blaze
    ${PROJECT_BINARY_DIR}/submodules/libxc
)
target_link_libraries(qcGPW_Harness PUBLIC qcUmbrella qchem qcMath qcCommon)

# THE PLANE-WAVE SCF HARNESS (qchem.Tests.PW_Harness): the prototype loop + the framework driver PW runs on
# until the facade grows a basis-family axis.  Over qcPW_TestFields (src/BasisSet/Lattice/tests).
add_library(qcPW_Harness)
target_sources(qcPW_Harness PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES PW/Harness.C)
target_compile_definitions(qcPW_Harness PUBLIC ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET qcPW_Harness PROPERTY CXX_STANDARD 20)
target_include_directories(qcPW_Harness PRIVATE
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    ${PROJECT_SOURCE_DIR}/submodules/blaze
    ${PROJECT_BINARY_DIR}/submodules/libxc
)
target_link_libraries(qcPW_Harness PUBLIC qcPW_TestFields qcUmbrella qchem qcMath qcCommon)

add_executable(ITMain gtestmain.C
    # THE PW SCF GRID (its unit-level tests left for src/ in phase 5)
    PW/Model.C PW/Si.C PW/NaF.C PW/CsI.C PW/Boxes.C
    A_HF_U.C A_HF_P.C A_HF_dfPin.C A_DHF.C A_DFT.C A_DFT_U.C A_PP.C L_PP.C M_HF_U.C M_DFT.C M_Sym.C
    # THE GPW SCF GRID (doc/TestSuitePlan.md §4): directory per basis family, file per material, ONE exe.
    GPW/Si.C GPW/Al.C GPW/Na.C GPW/NaF.C GPW/MnO.C GPW/Boxes.C
    M_Calculation.C
    M_Umbrella.C
    )
# The QchemTester/TestAtom/TestDiracAtom/TestMolecule scaffold is RETIRED: every molecular test drives
# qchem::Calculation and every atom/Dirac/PP test (and scfrun) drives qchem::AtomCalculation.  The shared
# RelativeError / RelativeHF/DFT/DHFError oracle helpers live in qchem.PeriodicTable (R2.20) -- the
# per-target TestUtils module file set is gone.

#-Wno-unknown-pragmas
# -Wno-cpp  Suppress warning messages emitted by #warning directives.

target_compile_definitions(ITMain PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
target_compile_definitions(ITMain PRIVATE ${BasisDataPath})
set_property(TARGET ITMain PROPERTY CXX_STANDARD 20)

target_include_directories(ITMain PRIVATE
    ${PROJECT_SOURCE_DIR}/submodules/BSplinebasis/include
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    ${PROJECT_SOURCE_DIR}/submodules/blaze
    ${PROJECT_SOURCE_DIR}/submodules/tabulate/include
    ${PROJECT_SOURCE_DIR}/submodules/googletest/googletest/include
    ${PROJECT_BINARY_DIR}/submodules/libxc
)
target_link_libraries(ITMain PRIVATE qcGPW_Harness qcPW_Harness qcUmbrella qchem qcMath qcCommon xc blas gtest)



# add_executable(DiracHydrogen DiracHydrogen.C)
# target_include_directories(DiracHydrogen PRIVATE ${PROJECT_SOURCE_DIR}/submodules/blaze)
# target_link_libraries(DiracHydrogen PUBLIC qchem lapack blas)
# set_property(TARGET DiracHydrogen PROPERTY CXX_STANDARD 20)

# setup testing
include(GoogleTest)
# DISCOVERY_EXTRA_ARGS --gtest_output= : force CMake's TEXT-listing discovery.  CMake >= 3.30 lists the tests
# through --gtest_output=json, and googletest escapes every non-ASCII BYTE of a name as \u00XX, so a test named
# with a Greek Gamma came back double-encoded ("Î" + U+0093) and its ctest filter matched NOTHING -- gtest ran
# zero tests, exited 0, and ctest printed "Passed 0.02 sec".  Every Γ-named SCF test (48 of them) was silently
# skipped by every sweep from the 2026-09-15 renaming until 2026-09-20.  A trailing empty --gtest_output=
# overrides the JSON request (gtest takes the last one), no file is written, and the module falls back to
# parsing the plain listing, which carries the bytes through intact.  Check: `ctest -N | grep Γ` must list them
# with the Gamma itself, never "Î".
gtest_discover_tests(ITMain DISCOVERY_EXTRA_ARGS --gtest_output=)
