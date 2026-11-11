# Tests for the qchem::Calculation facades

add_executable(UTCalculation ../../Common/tests/gtestmain.C
    RunDiagnostics.C
    Deck.C              # the input deck's data layer: typed options <-> JSON, --set, numbered revisions (D-ENV step 6)
    Materials.C
    ValenceBasisGen.C   # the valence-basis generator (from IntegrationTests 2026-09-15, TE phase 5)
    SolidState.C        # saved periodic SCF states: save / exact resume / warm start / refusal (CK-1)
    MoleculeSetup.C     # the molecular run setup: spin promotion, valence count, refusals (step 7)
    CalculationFacade.C # the front-end face: the ctor observer, levels, the spin density (lane G)
    MoleculeSession.C   # the molecular potential energy surface: reference == single point, carried step == fresh (step 8)
    MolecularRelaxation.C # REAL relaxations: the optimizer over a MoleculeSession, HF/DZVP, against PySCF minima (OpenWork §2)
    SolidSetup.C        # the periodic run setup: the frozen basis recipe, provenance without a deck (step 6)
    SolidSession.C      # the solid potential energy surface: reference == single point, carried step == fresh (step 8)
    DegenerateFill.C    # the degenerate-level fill policy beside GDM: the teaching `choose` row, the `warn` row (D-DEGEN-FILL follow-up)
    AuditFootprint.C    # the audit RSS policy: the byte count, the calibration against its record, the refusal doors (OpenWork (e))
    ForceFacade.C       # W and the per-term finite difference at fixed D, gated on Ven / kinetic by libcint (GeometryOptimizationPlan §6)
    WyckoffMap.C        # the level-1 solid map (free Wyckoff parameters) on periodic model surfaces (GeometryOptimizationPlan §1)
    GeometryOptimizer.C # the geometry optimizer on analytic model surfaces only (doc/GeometryOptimizationPlan.md §1, ruling §12.3)
    MeshAudit.C         # the molecular XC-mesh audit at a fixed D (GeometryOptimizationPlan §14.3)
    FacadeMatrix.C      # the facade PARITY MATRIX: {atom, molecule, solid} x every capability, expected failures named (V-FACADE-MATRIX)
)

target_compile_definitions(UTCalculation PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS} DECKS_DIR="${PROJECT_SOURCE_DIR}/decks")
set_property(TARGET UTCalculation PROPERTY CXX_STANDARD 20)

target_include_directories(UTCalculation PRIVATE
    ${PROJECT_SOURCE_DIR}/src
    ${QCHEM_BLAZE_DIR}
    ${QCHEM_JSON_INCLUDE}
    ${QCHEM_GTEST_INCLUDE}
)

target_link_libraries(UTCalculation PRIVATE qchem qcAD gtest)

qchem_discover_tests(UTCalculation)
