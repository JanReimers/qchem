
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

add_executable(ITMain gtestmain.C
    EigenSolverUT.C
    PW/PlaneWaveDFT.C
    A_HF_U.C A_HF_P.C A_HF_dfPin.C A_DHF.C A_DFT.C A_DFT_U.C A_PP.C L_PP.C GPW_UT.C M_HF_U.C M_DFT.C M_Sym.C
    # THE GPW SCF GRID (doc/TestSuitePlan.md §4): directory per basis family, file per material, ONE exe.
    GPW/Si.C GPW/Al.C GPW/Na.C GPW/NaF.C GPW/MnO.C GPW/Boxes.C
    M_Calculation.C
    RealComplexTermsUT.C
    ValenceBasisGen_UT.C
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
target_link_libraries(ITMain PRIVATE qcGPW_Harness qcUmbrella qchem qcMath qcCommon xc blas gtest)



# add_executable(DiracHydrogen DiracHydrogen.C)
# target_include_directories(DiracHydrogen PRIVATE ${PROJECT_SOURCE_DIR}/submodules/blaze)
# target_link_libraries(DiracHydrogen PUBLIC qchem lapack blas)
# set_property(TARGET DiracHydrogen PROPERTY CXX_STANDARD 20)

# setup testing
include(GoogleTest)
gtest_discover_tests(ITMain)
