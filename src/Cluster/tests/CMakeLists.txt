# Tests for qcCluster

add_executable(UTCluster ../../Common/tests/gtestmain.C
    ClusterTests.C
)

target_compile_options(UTCluster PRIVATE ${CXX_COMMON_COMPILE_OPTIONS})
target_compile_definitions(UTCluster PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET UTCluster PROPERTY CXX_STANDARD 20)

target_include_directories(UTCluster PRIVATE
    # ${PROJECT_SOURCE_DIR}/submodules/blaze
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    # ${PROJECT_SOURCE_DIR}/submodules/tabulate/include
)

target_link_libraries(UTCluster PRIVATE qcCluster)

gtest_discover_tests(UTCluster)
