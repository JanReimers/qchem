# Tests for qcLASolver

add_executable(UTLASolver ../../Common/tests/gtestmain.C
    LASolver.C
)

target_compile_options(UTLASolver PRIVATE ${CXX_COMMON_COMPILE_OPTIONS})
target_compile_definitions(UTLASolver PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET UTLASolver PROPERTY CXX_STANDARD 20)

target_include_directories(UTLASolver PRIVATE
    ${PROJECT_SOURCE_DIR}/submodules/googletest/googletest/include
)

target_link_libraries(UTLASolver PRIVATE qcLASolver lapack gtest)

gtest_discover_tests(UTLASolver)
