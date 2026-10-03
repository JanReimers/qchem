# qchem::Calculation -- the molecular SCF front-door facade.  Part of the top-level `qchem`
# aggregate library (it leans on the basis/Hamiltonian/accelerator/iterator pieces already linked
# there).  The Imp #includes nlohmann/json for the two json factories, so the json include lands here.
target_sources(qchem PRIVATE Imp/Calculation.C Imp/AtomCalculation.C Imp/SolidCalculation.C Imp/ValenceBasisGen.C Imp/Materials.C Imp/SolidState.C Imp/Deck.C)
# qchem.Materials reads the `species` (valence-count) field of src/Structure/Data/materials.json at run time -- INTERIM,
# D-STRUCTDATA step 2 moves it to its own pseudopotentials.json.  The structure fields are qcStructure's (qchem.StructureData).
target_compile_definitions(qchem PRIVATE ${StructureDataPath})
# all_my_modules' base dir was set by the first contributor (src/SCFIterator); register this
# directory as an additional base so the file set will accept Calculation.C from here.
target_sources(qchem PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES
    BASE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}
    FILES
    Calculation.C
    AtomCalculation.C
    SolidCalculation.C
    SolidState.C
    ValenceBasisGen.C
    Materials.C
    Deck.C
)
target_include_directories(qchem PRIVATE
    ${PROJECT_SOURCE_DIR}/src                       # forward.H (the unit-test friend declarations)
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
)

add_subdirectory(tests)
