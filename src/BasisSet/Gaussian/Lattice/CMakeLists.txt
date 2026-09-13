# src/BasisSet/Gaussian/Lattice — the G = T tier of the Gaussian engine: Bloch sums of the molecular
# Gaussians.  LatticeSum1E (the periodic-1E capability a PG basis realises), LatticeScreener, the spherical
# lattice view, and GPW (GPW_Evaluator + GPW_IBS: Gaussian orbitals over a plane-wave density fit -- the
# frozen family pairing of plan §5, hence qcGaussian_BS -> qcPlaneWave_BS).
# Modules: qchem.BasisSet.Gaussian.Lattice.*
target_sources(qcGaussian_BS PRIVATE
    Imp/SphericalLatticeView.C
    Imp/GPW_Evaluator.C
    Imp/GPW_IBS.C
)
target_sources(qcGaussian_BS PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES
    LatticeSum1E.C
    LatticeScreener.C
    SphericalLatticeView.C
    GPW_Evaluator.C
    GPW_IBS.C
)
target_include_directories(qcGaussian_BS PRIVATE ${PROJECT_SOURCE_DIR}/submodules/blaze)
