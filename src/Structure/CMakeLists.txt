add_library(qcStructure) 
target_sources(qcStructure PRIVATE 
    Imp/Atom.C 
    Imp/Structure.C 
    Imp/Molecule.C 
    Imp/UnitCell.C 
    Imp/Lattice.C
    Internal/Imp/AtomMesh.C 
    Internal/Imp/MoleculeMesh.C 
    )
target_sources(qcStructure PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES
    Structure.C
    UnitCell.C
    ReciprocalLattice.C
    KMesh.C
    Lattice.C
    Internal/UnitCellMesh.C
    Internal/AtomMesh.C
    Internal/MoleculeMesh.C
    )
target_include_directories(qcStructure PRIVATE 
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    ${PROJECT_SOURCE_DIR}/submodules/blaze
)
set_property(TARGET qcStructure PROPERTY CXX_STANDARD 20)
target_compile_options(qcStructure PRIVATE ${CXX_COMMON_COMPILE_OPTIONS})
target_compile_definitions(qcStructure PUBLIC ${CXX_COMMON_COMPILE_DEFINITIONS})
target_link_libraries(qcStructure PUBLIC qcMesh qcCommon gtest)

add_subdirectory(tests)