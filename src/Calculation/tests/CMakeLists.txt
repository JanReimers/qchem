# Tests for the qchem::Calculation facades

add_executable(UTCalculation ../../Common/tests/gtestmain.C
    RunDiagnostics.C
)

target_compile_definitions(UTCalculation PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET UTCalculation PROPERTY CXX_STANDARD 20)

target_include_directories(UTCalculation PRIVATE
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/submodules/blaze
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    ${PROJECT_SOURCE_DIR}/submodules/googletest/googletest/include
)

target_link_libraries(UTCalculation PRIVATE qchem gtest)

gtest_discover_tests(UTCalculation)
