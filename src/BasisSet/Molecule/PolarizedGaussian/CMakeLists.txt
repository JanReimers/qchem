target_sources(qcMolecule_BS PRIVATE
Imp/BasisSet.C
Imp/Symmetry.C
Imp/IrrepBasisSet.C
Imp/Fit_IBS.C
Imp/Orbital_IBS.C
Imp/SymmetryAdapt.C
Internal/Imp/Omega.C
Internal/Imp/Block.C
Internal/Imp/PGData.C
Internal/Imp/Polarization.C
MnD/Imp/Hermite2.C
MnD/Imp/Hermite1.C
Readers/Imp/Gaussian94.C
Radial/Imp/GaussianH3.C
Radial/Imp/GaussianRF.C
)
target_sources(qcMolecule_BS PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES
    Internal/PGData.C
    BasisSet.C
    Symmetry.C
    SymmetryAdapt.C
    Reader.C

    Internal/Omega.C
    Internal/Block.C
    Internal/Polarization.C
    MnD/Hermite1.C
    MnD/Hermite2.C
    MnD/Hermite3.C

    Radial/GData.C
    Radial/GaussianH3.C
    Radial/GaussianRF.C
    Readers/Gaussian94.C
)
target_include_directories(qcMolecule_BS PRIVATE ${PROJECT_SOURCE_DIR}/submodules/blaze)

 