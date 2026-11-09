# Tests for qcAD (the facade-level XC checks, which need a converged run, are in src/Calculation/tests/ForceFacade.C)
add_executable(UTAD ../../Common/tests/gtestmain.C
    BeckeAD.C
    XCKernelAD.C
)
target_compile_definitions(UTAD PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET UTAD PROPERTY CXX_STANDARD 20)
target_include_directories(UTAD PRIVATE
    ${PROJECT_SOURCE_DIR}/submodules/blaze
    ${PROJECT_SOURCE_DIR}/submodules/googletest/googletest/include
)
target_link_libraries(UTAD PRIVATE qcAD qcMesh gtest)
qchem_discover_tests(UTAD)
