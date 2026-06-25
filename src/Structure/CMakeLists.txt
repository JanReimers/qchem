add_library(qcStructure) 
target_sources(qcStructure PRIVATE 
    Imp/Atom.C 
    Imp/Structure.C 
    Imp/Molecule.C 
    Imp/UnitCell.C
    Lattice_3D/Imp/Lattice_3D.C
    Imp/MolecularMesh1.C
    )
target_sources(qcStructure PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES
    Structure.C
    UnitCell.C
    Ewald.C
    Lattice_3D/ReciprocalLattice.C
    Lattice_3D/KMesh.C
    Lattice_3D/Lattice_3D.C
    MolecularMesh1.C
    )
target_include_directories(qcStructure PRIVATE 
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    ${PROJECT_SOURCE_DIR}/submodules/blaze
)
set_property(TARGET qcStructure PROPERTY CXX_STANDARD 20)
target_compile_options(qcStructure PRIVATE ${CXX_COMMON_COMPILE_OPTIONS})
target_compile_definitions(qcStructure PUBLIC ${CXX_COMMON_COMPILE_DEFINITIONS})
target_link_libraries(qcStructure PUBLIC qcMesh qcMesh1 qcCommon)

add_subdirectory(tests)