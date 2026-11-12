# Tests for qcChargeDensity

add_executable(UTChargeDensity ../../Common/tests/gtestmain.C
   Version.C
   DensityMatrixDistance.C
   IonicCharges.C
   IonicSADSeed.C
   KerkerMix.C
   JointPulay.C
   AtomicDensity.C
   SinglesSampler.C
)

target_compile_definitions(UTChargeDensity PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET UTChargeDensity PROPERTY CXX_STANDARD 20)

target_include_directories(UTChargeDensity PRIVATE
    ${QCHEM_BLAZE_DIR}
    ${QCHEM_GTEST_INCLUDE}
)

# qcLattice_BS supplies the plane-wave block (SeedCD); qcChargeDensity the density types.
# lapack/blas resolve the linear-algebra symbols qcBasisSet/qcLASolver pull in transitively.
target_link_libraries(UTChargeDensity PRIVATE qcChargeDensity qcLattice_BS lapack blas gtest)

qchem_discover_tests(UTChargeDensity)
