target_sources(qcMolecule_BS PRIVATE
Imp/BasisSet.C
Imp/IrrepBasisSet.C
Imp/Orbital_IBS.C
Imp/Fit_IBS.C
Imp/Symmetry.C
)
target_sources(qcMolecule_BS PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES
    BasisSet.C
    Symmetry.C
)
target_include_directories(qcMolecule_BS PRIVATE ${PROJECT_SOURCE_DIR}/submodules/blaze)
