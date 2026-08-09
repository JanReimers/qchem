# Tests for qcSCFAccelerator

add_executable(UTSCFAccelerator ../../Common/tests/gtestmain.C
    SCFAccelerator.C
)

target_compile_options(UTSCFAccelerator PRIVATE ${CXX_COMMON_COMPILE_OPTIONS})
target_compile_definitions(UTSCFAccelerator PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET UTSCFAccelerator PROPERTY CXX_STANDARD 20)

target_include_directories(UTSCFAccelerator PRIVATE
    ${PROJECT_SOURCE_DIR}/submodules/googletest/googletest/include
    ${PROJECT_SOURCE_DIR}/submodules/blaze
)

target_link_libraries(UTSCFAccelerator PRIVATE qcSCFAccelerator lapack gtest)

gtest_discover_tests(UTSCFAccelerator)
