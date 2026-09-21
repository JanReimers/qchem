# Tests for qcHamiltonian

add_executable(UTHamiltonian ../../Common/tests/gtestmain.C
   Libxc.C
   LDA_XC.C
   EagerRefresh.C
   EnergyBreakdown.C
   RealComplexTerms.C   # the real/complex term twins on the GPW basis (from IntegrationTests 2026-09-15, TE phase 5)
   GPW_XC_FD.C          # the E_xc/H_xc FD-consistency probes (idem)
   Orthogonalize.C      # LASolver orthogonalization through a real atomic Hamiltonian (idem)
   PW_Terms.C           # the plane-wave term tests against the basis's direct-grid oracles (idem)
   Hubbard.C            # DFT+U: the LowdinProjector forward/adjoint pair (2026-09-19)
)

target_compile_definitions(UTHamiltonian PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET UTHamiltonian PROPERTY CXX_STANDARD 20)

target_include_directories(UTHamiltonian PRIVATE
    ${PROJECT_SOURCE_DIR}/submodules/blaze
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    ${PROJECT_SOURCE_DIR}/submodules/BSplinebasis/include
    ${PROJECT_SOURCE_DIR}/submodules/googletest/googletest/include
    ${PROJECT_BINARY_DIR}/submodules/libxc
)

# xc supplies the raw libxc C API exercised by Libxc.C and the oracle side of LDA_XC.C;
# qcHamiltonian supplies the production VWN_Correlation under test.
# qcLattice_BS + qcChargeDensity + qcFactory: the moved tests build GPW / periodic densities / a hydrogen atom
# basis as FIXTURES for the terms under test.  Siblings and leaves of qcHamiltonian, never anything above it.
target_link_libraries(UTHamiltonian PRIVATE qcPW_TestFields qcHamiltonian qcLattice_BS qcChargeDensity qcFactory xc lapack blas gtest)

qchem_discover_tests(UTHamiltonian)
