
link_directories( 
    $<$<CONFIG:Debug>:  submodules/wignerSymbols/build submodules/libxc/build>
    $<$<CONFIG:Release>:  submodules/wignerSymbols/build submodules/libxc/build>
)

add_executable(UTMain gtestmain.C QchemTesterImp.C 
    EigenSolverUT.C 
    A_HF_U.C A_HF_P.C A_DHF.C A_HF_U.C A_DFT.C A_DFT_U.C M_PG_U.C M_PG_Sym.C M_PG1_vs_PG.C
    )
target_sources(UTMain PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES QchemTester.C)

#-Wno-unknown-pragmas
# -Wno-cpp  Suppress warning messages emitted by #warning directives.

target_compile_options(UTMain PRIVATE ${CXX_COMMON_COMPILE_OPTIONS})
target_compile_definitions(UTMain PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
target_compile_definitions(UTMain PRIVATE ${BasisDataPath})
set_property(TARGET UTMain PROPERTY CXX_STANDARD 20)

target_link_options(UTMain PRIVATE -pg)
target_include_directories(UTMain PRIVATE
    ${PROJECT_SOURCE_DIR}/submodules/BSplinebasis/include
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    ${PROJECT_SOURCE_DIR}/submodules/blaze
    ${PROJECT_SOURCE_DIR}/submodules/tabulate/include
    ${PROJECT_BINARY_DIR}/submodules/libxc
    ${PROJECT_SOURCE_DIR}/submodules/wignerSymbols/include
)
target_link_libraries(UTMain PUBLIC qchem qcCommon xc)

# scfrun: command-line SCF driver for accelerator tuning (no recompile).
add_executable(scfrun scfrun.C QchemTesterImp.C)
target_sources(scfrun PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES QchemTester.C )
target_compile_options(scfrun PRIVATE ${CXX_COMMON_COMPILE_OPTIONS})
target_compile_definitions(scfrun PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET scfrun PROPERTY CXX_STANDARD 20)
target_include_directories(scfrun PRIVATE
    ${PROJECT_SOURCE_DIR}/submodules/BSplinebasis/include
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    ${PROJECT_SOURCE_DIR}/submodules/blaze
    ${PROJECT_SOURCE_DIR}/submodules/tabulate/include
)
target_link_libraries(scfrun PUBLIC qchem qcCommon)


# add_executable(DiracHydrogen DiracHydrogen.C)
# target_include_directories(DiracHydrogen PRIVATE ${PROJECT_SOURCE_DIR}/submodules/blaze)
# target_link_libraries(DiracHydrogen PUBLIC qchem lapack blas)
# set_property(TARGET DiracHydrogen PROPERTY CXX_STANDARD 20)

# setup testing
enable_testing()
include(CTest)
include(GoogleTest)
# gtest_discover_tests(UTMain)
