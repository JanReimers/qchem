# Tests for qcChargeDensity

add_executable(UTChargeDensity ../../Common/tests/gtestmain.C
   Version.C
   IonicCharges.C
)

target_compile_options(UTChargeDensity PRIVATE ${CXX_COMMON_COMPILE_OPTIONS})
target_compile_definitions(UTChargeDensity PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET UTChargeDensity PROPERTY CXX_STANDARD 20)

target_include_directories(UTChargeDensity PRIVATE
    ${PROJECT_SOURCE_DIR}/submodules/blaze
)

# qcLattice_BS supplies the plane-wave block (SeedCD); qcChargeDensity the density types.
# lapack/blas resolve the linear-algebra symbols qcBasisSet/qcLASolver pull in transitively.
target_link_libraries(UTChargeDensity PRIVATE qcChargeDensity qcLattice_BS lapack blas gtest)

gtest_discover_tests(UTChargeDensity)
