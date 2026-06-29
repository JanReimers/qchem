# Tests for qcPseudopotential

add_executable(UTPseudopotential ../../Common/tests/gtestmain.C
   SeparableViews.C
)

target_compile_options(UTPseudopotential PRIVATE ${CXX_COMMON_COMPILE_OPTIONS})
target_compile_definitions(UTPseudopotential PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET UTPseudopotential PROPERTY CXX_STANDARD 20)

target_include_directories(UTPseudopotential PRIVATE
    ${PROJECT_SOURCE_DIR}/submodules/blaze
)

# qcMesh supplies the radial mesh for the Bessel-transform quadrature.
target_link_libraries(UTPseudopotential PRIVATE qcPseudopotential qcMesh gtest)

gtest_discover_tests(UTPseudopotential)
