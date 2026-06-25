# Tests for qcLattice_BS (true plane-wave / APW / LAPW UNIT tests -- no SCF stack).
#
# This app links ONLY qcLattice_BS and its down-DAG closure (qcBasisSet, qcStructure, qcSymmetry,
# qcLASolver, qcMesh, qcCommon).  It is deliberately NOT allowed to link anything higher in the DAG
# (Hamiltonian / ChargeDensity / Orbitals / WaveFunction / SCFIterator) -- so the linker structurally
# guarantees these stay true unit tests.  The self-consistent (SCF) integration tests that DO need the
# upper stack live in UnitTests/PlaneWaveDFTUT.C instead.

add_executable(UTLattice_3D_BS ../../Atom/tests/gtestmain.C
    PlaneWaveUT.C
    BandStructureUT.C
    APW_UT.C
    LAPW_UT.C
)

target_compile_options(UTLattice_3D_BS PRIVATE ${CXX_COMMON_COMPILE_OPTIONS})
target_compile_definitions(UTLattice_3D_BS PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS} ${LatticeDataPath})
set_property(TARGET UTLattice_3D_BS PROPERTY CXX_STANDARD 20)

target_include_directories(UTLattice_3D_BS PRIVATE
    ${PROJECT_SOURCE_DIR}/submodules/blaze
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
)

target_link_libraries(UTLattice_3D_BS PRIVATE qcLattice_BS lapack gtest)

gtest_discover_tests(UTLattice_3D_BS)
