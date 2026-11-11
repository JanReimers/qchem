# Tests for qcOptimize: the steppers on ANALYTIC model costs (the geometry-specific gates -- maps, symmetry, the PES -- stay in UTCalculation).
add_executable(UTOptimize ../../Common/tests/gtestmain.C
    Minimize.C
)
target_compile_definitions(UTOptimize PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET UTOptimize PROPERTY CXX_STANDARD 20)
target_include_directories(UTOptimize PRIVATE
    ${QCHEM_BLAZE_DIR}
    ${QCHEM_GTEST_INCLUDE}
)
target_link_libraries(UTOptimize PRIVATE qcOptimize gtest)
qchem_discover_tests(UTOptimize)
