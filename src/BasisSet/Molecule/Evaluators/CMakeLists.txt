# Molecular evaluators (Goal B: lift the atom Evaluator pattern to PG).  Copied from
# Atom/Evaluators and diverging; interface-only for now (no Imp yet).
target_sources(qcMolecule_BS PRIVATE
    Internal/Boys/Imp/AuxillaryFJ.C
    Internal/Boys/Imp/AuxillaryData.C
    Internal/MnD/Imp/Triangle3D.C
    Internal/MnD/Imp/RNLM.C
)
target_sources(qcMolecule_BS PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES
    Evaluator.C
    PG_Cart_MnD/Evaluator.C

    Internal/Boys/AuxillaryFJ.C
    Internal/MnD/Index3.C
    Internal/MnD/Triangle3D.C
    Internal/MnD/RNLM.C
    Internal/MnD/MnD.C
)
target_include_directories(qcMolecule_BS PRIVATE ${PROJECT_SOURCE_DIR}/submodules/blaze)
