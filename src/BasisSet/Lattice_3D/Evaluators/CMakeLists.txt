# src/BasisSet/Lattice_3D/Evaluators — the GPW orbital evaluator (Bloch sums of Gaussians over a PW
# density grid).  The plane-wave engine it builds on lives in ../../PlaneWave since plan step 1a; this
# one leaves for the Gaussian library at step 1b.
target_sources(qcLattice_BS PRIVATE
    GPW/Imp/Evaluator.C
)
target_sources(qcLattice_BS PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES
    GPW/Evaluator.C
)
