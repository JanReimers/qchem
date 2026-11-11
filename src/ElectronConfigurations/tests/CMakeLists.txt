# Tests for qcElConfig

add_executable(UTElConfig ../../Common/tests/gtestmain.C
   ElectronConfigurations.C
   ElectronConfigurations_Dirac.C
   OccupationState.C
   ResponseWeight.C    # the occupancy rules' response weights (LinearResponsePlan E1)
)

target_compile_definitions(UTElConfig PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET UTElConfig PROPERTY CXX_STANDARD 20)

target_include_directories(UTElConfig PRIVATE
    # ${QCHEM_BLAZE_DIR}
    # ${QCHEM_JSON_INCLUDE}
    ${QCHEM_TABULATE_INCLUDE}
    ${QCHEM_GTEST_INCLUDE}
)

target_link_libraries(UTElConfig PRIVATE qcElConfig gtest)

qchem_discover_tests(UTElConfig)
