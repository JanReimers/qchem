# src/BasisSet/Lattice_3D/Evaluators — plane-wave grid evaluators (the lattice analog of
# Molecule/Evaluators).  One evaluator today (PW); expected to grow (APW/LAPW/GPW ...).
target_sources(qcLattice_BS PRIVATE
    PW/Imp/Evaluator.C
)
target_sources(qcLattice_BS PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES
    PW/Evaluator.C
)
