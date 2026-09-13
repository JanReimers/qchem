# src/BasisSet/Gaussian/Point — the G = P (point-group) tier of the Gaussian engine: the molecular IBS
# mixins, PG_Cart / PG_Spherical / PG_LibCint, the SALC container, the factory and the basis-file readers.
# Modules: qchem.BasisSet.Gaussian.Point.*  (must not import a lattice -- plan step 3 audit).
target_sources(qcGaussian_BS PRIVATE
    Imp/Factory.C
    Imp/SymmetryAdaptedBasisSet.C
    Imp/BasisFiles.C
    Readers/Imp/Gaussian94.C
)
target_sources(qcGaussian_BS PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES
    Factory.C
    SymmetryAdaptedBasisSet.C
    IrrepBasisSet.C
    BasisFiles.C
    Reader.C
    Readers/Gaussian94.C
)
add_subdirectory(PG_Cart)
add_subdirectory(PG_Spherical)
add_subdirectory(PG_LibCint)
