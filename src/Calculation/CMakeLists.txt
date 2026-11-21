# The SCF facades: qchem::Calculation (the abstract base), qchem::MoleculeCalculation (the molecular front door),
# qchem::AtomCalculation, qchem::SolidCalculation.  Part of the top-level `qchem`
# aggregate library (it leans on the basis/Hamiltonian/accelerator/iterator pieces already linked
# there).  The Imp #includes nlohmann/json for the two json factories, so the json include lands here.
target_sources(qchem PRIVATE Imp/MoleculeCalculation.C Imp/Calculation.C
    Imp/AtomCalculation.C Imp/SolidCalculation.C Imp/SolidSetup.C Imp/MoleculeSetup.C Imp/AtomSetup.C Imp/ValenceBasisGen.C Imp/Materials.C Imp/SolidState.C Imp/SCFState.C Imp/StateView.C Imp/Deck.C Imp/PseudoAtomOrbital.C Imp/AuditFootprint.C Imp/SetupDecisions.C Imp/BasisMenu.C)
# qchem.Materials reads the `pp` (pseudopotential set name) field of src/Structure/Data/materials.json at run time; the set is qchem.Pseudopotential.Set's.
# The structure fields are qcStructure's (qchem.StructureData).
target_compile_definitions(qchem PRIVATE ${StructureDataPath})
# all_my_modules' base dir was set by the first contributor (src/SCFIterator); register this
# directory as an additional base so the file set will accept MoleculeCalculation.C from here.
target_sources(qchem PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES
    BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
    FILES
    Calculation.C
    MoleculeCalculation.C
    AtomCalculation.C
    SolidCalculation.C
    SolidSetup.C
    MoleculeSetup.C
    PotentialEnergySurface.C
    GeometryOptimizer.C
    SCFState.C
    StateView.C
    WyckoffMap.C
    SolidState.C
    ValenceBasisGen.C
    PseudoAtomOrbital.C
    Materials.C
    Deck.C
    SCFFailure.C
    SetupRefusal.C
    Setup.C                               # the common run setup (V-FACADE-PARITY P1)
    AtomSetup.C                           # the radial atom's run setup (V-FACADE-PARITY P1, R3)
    SCFRecipeReport.C
    SCFRecipe.C
    DegenerateFill.C
    Policy.C
    GridAudit.C
    BasisCoverage.C
    BasisMenu.C                           # which orbital bases a structure admits: the front ends' menu (lane G 2026-11-21)
    XCRhoPolicy.C
    SeedPolicy.C
    RunOptions.C                          # THE one options tree (V-FACADE-PARITY P0b-2)
    SetupDecisions.C                      # the structure-neutral setup decisions (V-FACADE-SETUP)
    AuditFootprint.C
    PhiTablesPolicy.C                     # the δ fit's XC orbital tables: dense or shared images, a memory-vs-speed choice (lane S 2026-11-19)
    ERI4MemoryPolicy.C                    # the ERI4 RAM cap as a setup policy (IntegralCacheOwnershipPlan S6, ruling 6)
)
target_include_directories(qchem PRIVATE
    ${PROJECT_SOURCE_DIR}/src                       # forward.H (the unit-test friend declarations)
    ${QCHEM_JSON_INCLUDE}
)

add_subdirectory(tests)
