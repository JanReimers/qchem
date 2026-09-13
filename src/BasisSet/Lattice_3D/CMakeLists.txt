# src/BasisSet/Lattice_3D — the G=T basis-set layer above the engines (doc/BasisSetTaxonomyPlan.md §3).
# Plan step 1a (2026-09-13) split the plane-wave ENGINE out into ../PlaneWave (qcPlaneWave_BS); what is
# here now is the GPW engine + GPW_IBS (leaving for the Gaussian library at step 1b) and the T container
# tier (BasisSet / BandStructure / APW / LAPW) that stays.  See doc/OldPlans/PlaneWavePlan.md for the history.
add_library(qcLattice_BS)
target_sources(qcLattice_BS PRIVATE
    Imp/GPW_IBS.C
    Imp/APW_IBS.C
    Imp/LAPW_IBS.C
    Imp/BasisSet.C
)
target_sources(qcLattice_BS PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES
    GPW_IBS.C
    APW_IBS.C
    LAPW_IBS.C
    BandStructure.C
    BasisSet.C
)
target_include_directories(qcLattice_BS PRIVATE
    ${PROJECT_SOURCE_DIR}/submodules/blaze
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    ${PROJECT_SOURCE_DIR}/src
)
set_property(TARGET qcLattice_BS PROPERTY CXX_STANDARD 20)
target_compile_definitions(qcLattice_BS PUBLIC ${CXX_COMMON_COMPILE_DEFINITIONS})
# qcPlaneWave_BS: the plane-wave engine (PlaneWave_IBS, the PW fit basis, the FFT/Poisson grid) -- split
# out at plan step 1a; PUBLIC so every client that reached those modules through qcLattice_BS still does.
# qcMolecule_BS: GPW's periodic Gaussian 1E integrals delegate to the molecular Gaussian basis's
# Molecule::LatticeSum1E capability.  No cycle (qcMolecule_BS never references Lattice).
# V1.2 (2026-09-13): qcPseudopotential is NO LONGER a dependency -- the species-field faces the PP models implement
# live in qcBasisSet (Orbital_PP_IBS / SpeciesField), and the edge now runs the other way.
target_link_libraries(qcLattice_BS PUBLIC qcBasisSet qcStructure qcSymmetry qcLASolver qcPlaneWave_BS qcMolecule_BS)

add_subdirectory(Evaluators)
add_subdirectory(tests)
