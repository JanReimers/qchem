# Tests for qcCommon

add_executable(qcCommon_tests gtestmain.C
    PeriodicTable.C
)

target_compile_options(qcCommon_tests PRIVATE ${CXX_COMMON_COMPILE_OPTIONS})
target_compile_definitions(qcCommon_tests PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET qcCommon_tests PROPERTY CXX_STANDARD 20)

target_include_directories(qcCommon_tests PRIVATE
    ${PROJECT_SOURCE_DIR}/submodules/blaze
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    ${PROJECT_SOURCE_DIR}/submodules/tabulate/include
)

target_link_libraries(qcCommon_tests PRIVATE qcCommon gtest)

gtest_discover_tests(qcCommon_tests)
