# Tests for qcCommon

add_executable(UTCommon gtestmain.C
    PeriodicTable.C
)

target_compile_options(UTCommon PRIVATE ${CXX_COMMON_COMPILE_OPTIONS})
target_compile_definitions(UTCommon PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET UTCommon PROPERTY CXX_STANDARD 20)

target_include_directories(UTCommon PRIVATE
    # ${PROJECT_SOURCE_DIR}/submodules/blaze
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    # ${PROJECT_SOURCE_DIR}/submodules/tabulate/include
)

target_link_libraries(UTCommon PRIVATE qcCommon gtest)

gtest_discover_tests(UTCommon)
