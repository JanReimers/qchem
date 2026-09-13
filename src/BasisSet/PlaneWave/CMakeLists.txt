# src/BasisSet/PlaneWave — the PLANE-WAVE engine library (doc/BasisSetTaxonomyPlan.md §3, split out of
# qcLattice_BS at step 1a, 2026-09-13).  FAMILY = exponential e^{i(k+G).r}; G = T (the Bloch label).
#
# Contents: the PW_Evaluator / PW_Grid_Evaluator engines + the shared FFT/Poisson PeriodicGridEvaluator,
# the orbital PlaneWave_IBS and the auxiliary PlaneWaveFit_IBS (density + Vxc fit faces), and the
# Internal G-vector / k+G helpers.  The G=T SPEC these engines satisfy (isLattice_{1E,DFT}_Evaluator +
# the Lattice::*_IBS mixins) lives in the core (qchem.BasisSet.Lattice_IBS), not here.
#
# A LEAF over the core: links qcBasisSet + qcStructure + qcSymmetry only -- NOTHING from qcMolecule_BS or
# the qcLattice_BS container tier (measured 2026-09-13; the GPW engine depends on THIS library, never the
# reverse).  Module names keep their qchem.BasisSet.Lattice_3D.* spelling until plan step 2 renames them
# onto the two axes (qchem.BasisSet.PlaneWave.*).
add_library(qcPlaneWave_BS)
target_sources(qcPlaneWave_BS PRIVATE
    Evaluators/Imp/PeriodicGridEvaluator.C
    Evaluators/Imp/Evaluator.C
    Imp/PlaneWave_IBS.C
)
target_sources(qcPlaneWave_BS PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES
    Internal/GVectors.C
    Internal/KPlusG.C
    Evaluators/PeriodicGridEvaluator.C
    Evaluators/Evaluator.C
    PlaneWave_IBS.C
    PlaneWaveFit_IBS.C
)
target_include_directories(qcPlaneWave_BS PRIVATE
    ${PROJECT_SOURCE_DIR}/submodules/blaze
    ${PROJECT_SOURCE_DIR}/src
)
set_property(TARGET qcPlaneWave_BS PROPERTY CXX_STANDARD 20)
target_compile_definitions(qcPlaneWave_BS PUBLIC ${CXX_COMMON_COMPILE_DEFINITIONS})
target_link_libraries(qcPlaneWave_BS PUBLIC qcBasisSet qcStructure qcSymmetry)
