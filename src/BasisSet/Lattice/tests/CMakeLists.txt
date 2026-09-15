# Tests for qcLattice_BS (true plane-wave / APW / LAPW UNIT tests -- no SCF stack).
#
# This app links ONLY qcLattice_BS and its down-DAG closure (qcBasisSet, qcStructure, qcSymmetry,
# qcLASolver, qcMesh, qcCommon).  It is deliberately NOT allowed to link anything higher in the DAG
# (Hamiltonian / ChargeDensity / Orbitals / WaveFunction / SCFIterator) -- so the linker structurally
# guarantees these stay true unit tests.  The self-consistent (SCF) integration tests that DO need the
# upper stack live in IntegrationTests/{GPW,PW}/ instead; the two GPW XC finite-difference probes that need
# the XC functional objects live in src/Hamiltonian/tests/GPW_XC_FD.C for the same reason.

# THE PLANE-WAVE TEST FIXTURES (qchem.Tests.PW_Fields): direct-grid oracles + PWFixture, Hamiltonian-free.  A
# library so the PW term tests (UTHamiltonian) and the PW SCF grid (IntegrationTests/PW) share one copy.
add_library(qcPW_TestFields)
target_sources(qcPW_TestFields PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES PW_Fields.C)
target_compile_definitions(qcPW_TestFields PUBLIC ${CXX_COMMON_COMPILE_DEFINITIONS})
set_property(TARGET qcPW_TestFields PROPERTY CXX_STANDARD 20)
target_include_directories(qcPW_TestFields PRIVATE ${PROJECT_SOURCE_DIR}/submodules/blaze ${PROJECT_SOURCE_DIR}/submodules/json/single_include)
target_link_libraries(qcPW_TestFields PUBLIC qcLattice_BS qcPseudopotential qcFitting)

add_executable(UTLattice_BS ../../Radial/tests/gtestmain.C
    PlaneWaveUT.C
    BandStructureUT.C
    APW_UT.C
    LAPW_UT.C
    GTH_UT.C
    GMapUT.C
    ScreenedIntegratorUT.C
    RealComplexUT.C
    GPW_UT.C            # the GPW basis/evaluator unit tests (moved from IntegrationTests 2026-09-15, TE phase 5)
    PlaneWave_FieldsUT.C  # the PW basis / FFT / Poisson unit tests (idem)
)

target_compile_definitions(UTLattice_BS PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS})  # GTH data path now baked into qcPseudopotential
set_property(TARGET UTLattice_BS PROPERTY CXX_STANDARD 20)

target_include_directories(UTLattice_BS PRIVATE
    ${PROJECT_SOURCE_DIR}/submodules/blaze
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    ${PROJECT_SOURCE_DIR}/submodules/googletest/googletest/include
)

target_link_libraries(UTLattice_BS PRIVATE qcPW_TestFields qcLattice_BS qcPseudopotential qcFitting lapack gtest)   # the tests use the GTH models directly; qcFitting = OrthogonalFit (GPW_UT's per-representation fit gates), a leaf

gtest_discover_tests(UTLattice_BS)
