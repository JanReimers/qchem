# Tests for qcLASolver

add_executable(UTLASolver ../../Common/tests/gtestmain.C
    LASolver.C
    OverlapConditioning.C
    Krylov.C
)

target_compile_definitions(UTLASolver PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET UTLASolver PROPERTY CXX_STANDARD 20)

target_include_directories(UTLASolver PRIVATE
    ${QCHEM_GTEST_INCLUDE}
)

target_link_libraries(UTLASolver PRIVATE qcLASolver lapack gtest)

qchem_discover_tests(UTLASolver)
