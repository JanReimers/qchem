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
    SolidSetup.C        # the periodic run setup: the frozen basis recipe, provenance without a deck (step 6)
    SolidSession.C      # the solid potential energy surface: reference == single point, carried step == fresh (step 8)
    AuditFootprint.C    # the audit RSS policy: the byte count, the calibration against its record, the refusal doors (OpenWork (e))
    ForceFacade.C       # W and the per-term finite difference at fixed D, gated on Ven / kinetic by libcint (GeometryOptimizationPlan §6)
    GeometryOptimizer.C # the geometry optimizer on analytic model surfaces only (doc/GeometryOptimizationPlan.md §1, ruling §12.3)
    FacadeMatrix.C      # the facade PARITY MATRIX: {atom, molecule, solid} x every capability, expected failures named (V-FACADE-MATRIX)
)

target_compile_definitions(UTCalculation PRIVATE ${CXX_COMMON_COMPILE_DEFINITIONS} DECKS_DIR="${PROJECT_SOURCE_DIR}/decks")
set_property(TARGET UTCalculation PROPERTY CXX_STANDARD 20)

target_include_directories(UTCalculation PRIVATE
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/submodules/blaze
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    ${PROJECT_SOURCE_DIR}/submodules/googletest/googletest/include
)

target_link_libraries(UTCalculation PRIVATE qchem gtest)

qchem_discover_tests(UTCalculation)
