# Tests for qcPseudopotential

add_executable(UTPseudopotential ../../Common/tests/gtestmain.C
   SeparableViews.C
   PseudopotentialSet.C
)

target_compile_definitions(UTPseudopotential PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET UTPseudopotential PROPERTY CXX_STANDARD 20)

target_include_directories(UTPseudopotential PRIVATE
    ${QCHEM_BLAZE_DIR}
    ${QCHEM_GTEST_INCLUDE}
)

# qcMesh supplies the radial mesh for the Bessel-transform quadrature.
target_link_libraries(UTPseudopotential PRIVATE qcPseudopotential qcMesh gtest)

qchem_discover_tests(UTPseudopotential)
