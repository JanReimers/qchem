# src/BasisSet/Gaussian/PG_Cart — the raw Cartesian-Gaussian AO block: the family's G=1 SEED.  It answers
# the two seed questions the two constructions ask from OUTSIDE it -- GetAoShells (the P induction, done by
# Point/SymmetryAdaptedBasisSet) and LatticeSum1E (the T induction, done by Lattice/GPW) -- and forwards the
# lattice kernels to its PG_Cart_MnD engine base (where UnitCell is legitimate, plan §1.6 rule 1).  Untagged
# (qchem.BasisSet.Gaussian.PG_Cart.*) by the step-3 ruling; the §5 sequel (Point/Lattice specs + one thin
# class per (G, engine)) is where its dual role splits for real.
target_sources(qcGaussian_BS PRIVATE
Imp/BasisSet.C
Imp/Symmetry.C
Imp/IrrepBasisSet.C
Imp/Fit_IBS.C
Imp/Orbital_IBS.C
Imp/SymmetryAdapt.C
)
target_sources(qcGaussian_BS PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES
    BasisSet.C
    Symmetry.C
    SymmetryAdapt.C
)
target_include_directories(qcGaussian_BS PRIVATE ${PROJECT_SOURCE_DIR}/submodules/blaze)
