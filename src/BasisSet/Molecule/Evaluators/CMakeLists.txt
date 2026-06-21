# Molecular evaluators (Goal B: lift the atom Evaluator pattern to PG).  Copied from
# Atom/Evaluators and diverging; interface-only for now (no Imp yet).
target_sources(qcMolecule_BS PRIVATE
    Internal/Boys/Imp/AuxillaryFJ.C
    Internal/Boys/Imp/AuxillaryData.C
    Internal/MnD/Imp/Triangle3D.C
    Internal/MnD/Imp/RNLM.C
    # PG_Cart_MnD Cartesian glue (Imp)
    PG_Cart_MnD/Imp/PGData.C
    PG_Cart_MnD/Imp/Polarization.C
    PG_Cart_MnD/Imp/GaussianRF.C
    PG_Cart_MnD/Internal/Imp/Block.C
    PG_Cart_MnD/Internal/Imp/GaussianH3.C
    PG_Cart_MnD/Internal/Imp/Hermite1.C
    PG_Cart_MnD/Internal/Imp/Hermite2.C
)
target_sources(qcMolecule_BS PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES
    Evaluator.C
    PG_Cart_MnD/Evaluator.C

    Internal/Boys/AuxillaryFJ.C
    Internal/MnD/Index3.C
    Internal/MnD/Triangle3D.C
    Internal/MnD/RNLM.C
    Internal/MnD/MnD.C

    # PG_Cart_MnD Cartesian glue (module interfaces)
    PG_Cart_MnD/PGData.C
    PG_Cart_MnD/Polarization.C
    PG_Cart_MnD/GaussianRF.C
    PG_Cart_MnD/Internal/Block.C
    PG_Cart_MnD/Internal/Hermite.C     # Hermite1/2/3 + GaussianH3 (Imp split below)
)
target_include_directories(qcMolecule_BS PRIVATE ${PROJECT_SOURCE_DIR}/submodules/blaze)
