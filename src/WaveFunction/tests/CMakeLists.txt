# Tests for qcWaveFunction -- the level-table ROW builders (D10) and the R2.4 reservoir-charge check.
add_executable(UTWaveFunction ../../Common/tests/gtestmain.C
    EigenRowsUT.C
    ReservoirChargeUT.C
    DegenerateFillUT.C  # the kT = 0 fill kernel (D-DEGEN-FILL)
)
target_compile_definitions(UTWaveFunction PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET UTWaveFunction PROPERTY CXX_STANDARD 20)
target_include_directories(UTWaveFunction PRIVATE
    ${PROJECT_SOURCE_DIR}/submodules/googletest/googletest/include
)
target_link_libraries(UTWaveFunction PRIVATE qcWaveFunction qcOrbitals qcSymmetry gtest)
qchem_discover_tests(UTWaveFunction)
