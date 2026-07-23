# Tests for qcHamiltonian

add_executable(UTHamiltonian ../../Common/tests/gtestmain.C
   Libxc.C
   LDA_XC.C
)

target_compile_options(UTHamiltonian PRIVATE ${CXX_COMMON_COMPILE_OPTIONS})
target_compile_definitions(UTHamiltonian PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET UTHamiltonian PROPERTY CXX_STANDARD 20)

target_include_directories(UTHamiltonian PRIVATE
    ${PROJECT_SOURCE_DIR}/submodules/blaze
    ${PROJECT_SOURCE_DIR}/submodules/googletest/googletest/include
    ${PROJECT_BINARY_DIR}/submodules/libxc
)

# xc supplies the raw libxc C API exercised by Libxc.C and the oracle side of LDA_XC.C;
# qcHamiltonian supplies the production VWN_Correlation under test.
target_link_libraries(UTHamiltonian PRIVATE qcHamiltonian xc lapack blas gtest)

gtest_discover_tests(UTHamiltonian)
