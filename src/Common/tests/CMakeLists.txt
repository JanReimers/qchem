# Tests for qcCommon

add_executable(UTCommon gtestmain.C
    PeriodicTable.C
    Blaze.C
    FFT.C
    NUFFT.C
    Reporting.C
    Outcome.C
    RunPolicy.C
    Diagnostics.C
    HDF5.C   # the saved-SCF-state carrier (CK-1)
    PositionJet.C   # qchem.PositionJet (qcMath): the first-order jet in the atom positions
)

target_compile_definitions(UTCommon PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET UTCommon PROPERTY CXX_STANDARD 20)

target_include_directories(UTCommon PRIVATE
    # ${PROJECT_SOURCE_DIR}/submodules/blaze
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    ${PROJECT_SOURCE_DIR}/submodules/googletest/googletest/include
    # ${PROJECT_SOURCE_DIR}/submodules/tabulate/include
)

target_link_libraries(UTCommon PRIVATE qcMath qcCommon gtest)

qchem_discover_tests(UTCommon)
