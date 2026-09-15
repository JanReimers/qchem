# Command-line applications (not tests).

# scfrun: command-line SCF driver for accelerator tuning (no recompile).  Drives qchem::AtomCalculation
# directly (no QchemTester); the RelativeHF/DHFError oracle helpers come from qchem.PeriodicTable (R2.20)
# -- no test module imported.
add_executable(scfrun scfrun.C)
target_compile_definitions(scfrun PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET scfrun PROPERTY CXX_STANDARD 20)
target_include_directories(scfrun PRIVATE
    ${PROJECT_SOURCE_DIR}/submodules/BSplinebasis/include
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    ${PROJECT_SOURCE_DIR}/submodules/blaze
    ${PROJECT_SOURCE_DIR}/submodules/tabulate/include
)
target_link_libraries(scfrun PUBLIC qchem qcMath qcCommon blas)

# valgen: command-line valence-basis generator (GPWPlan1 step 1).  A thin arg-parser over the
# qchem.ValenceBasisGen library; dogfoods qchem.Reporting for its console output.
add_executable(valgen valgen.C)
target_compile_definitions(valgen PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET valgen PROPERTY CXX_STANDARD 20)
target_include_directories(valgen PRIVATE
    ${PROJECT_SOURCE_DIR}/submodules/BSplinebasis/include
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    ${PROJECT_SOURCE_DIR}/submodules/blaze
    ${PROJECT_SOURCE_DIR}/submodules/tabulate/include
)
target_link_libraries(valgen PUBLIC qchem qcMath qcCommon blas)

# gpwprobe: the GPW campaign INSTRUMENTS (supercell ladder, k sweep, NaF GDM-smear probe, the V2.6 Becke
# ladders, the MnO AFM-II campaign run) -- promoted out of IntegrationTests/GPW_SCF_UT.C's DISABLED_ tests
# (doc/TestSuitePlan.md §8).  Shares the gates' harness (qcGPW_Harness) so a probe and its gate build the
# same cell the same way.
add_executable(gpwprobe gpwprobe.C)
target_compile_definitions(gpwprobe PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET gpwprobe PROPERTY CXX_STANDARD 20)
target_include_directories(gpwprobe PRIVATE
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    ${PROJECT_SOURCE_DIR}/submodules/blaze
)
target_link_libraries(gpwprobe PUBLIC qcGPW_Harness qchem qcMath qcCommon blas)
