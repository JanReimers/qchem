# src/BasisSet/Lattice_3D — the G=T basis-set layer above the engines (doc/BasisSetTaxonomyPlan.md §3).
# Plan steps 1a/1b (2026-09-13) moved the two ENGINES out -- plane waves to ../PlaneWave (qcPlaneWave_BS),
# GPW to ../Molecule/Lattice (qcGaussian_BS).  What is here is the T CONTAINER tier: the ⊕_k BasisSet +
# its factories, BandStructure, and the APW/LAPW composites -- the one G whose container must see two
# engine libraries (§1.6 corollary).  See doc/OldPlans/PlaneWavePlan.md for the history.
add_library(qcLattice_BS)
target_sources(qcLattice_BS PRIVATE
    Imp/APW_IBS.C
    Imp/LAPW_IBS.C
    Imp/BasisSet.C
)
target_sources(qcLattice_BS PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES
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
# The two engines, both PUBLIC so every client that reached their modules through qcLattice_BS still does:
# qcPlaneWave_BS (PlaneWave_IBS, the PW fit basis, the FFT/Poisson grid; step 1a) and qcGaussian_BS
# (GPW_IBS + GPW_Evaluator since step 1b).  No cycle: neither engine references this container tier.
# V1.2 (2026-09-13): qcPseudopotential is NO LONGER a dependency -- the species-field faces the PP models implement
# live in qcBasisSet (Orbital_PP_IBS / SpeciesField), and the edge now runs the other way.
target_link_libraries(qcLattice_BS PUBLIC qcBasisSet qcStructure qcSymmetry qcLASolver qcPlaneWave_BS qcGaussian_BS)

add_subdirectory(tests)
