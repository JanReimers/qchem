add_library(qcCluster) 
target_sources(qcCluster PRIVATE 
    Imp/Atom.C 
    Imp/Cluster.C 
    Imp/Molecule.C 
    Imp/UnitCell.C 
    Imp/Lattice.C
    Internal/Imp/AtomMesh.C 
    Internal/Imp/MoleculeMesh.C 
    )
target_sources(qcCluster PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
    Cluster.C 
    UnitCell.C 
    Lattice.C 
    Internal/UnitCellMesh.C 
    Internal/AtomMesh.C 
    Internal/MoleculeMesh.C
    )
target_include_directories(qcCluster PRIVATE 
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
    ${PROJECT_SOURCE_DIR}/submodules/blaze
)
set_property(TARGET qcCluster PROPERTY CXX_STANDARD 20)
target_compile_options(qcCluster PRIVATE ${CXX_COMMON_COMPILE_OPTIONS})
target_compile_definitions(qcCluster PUBLIC ${CXX_COMMON_COMPILE_DEFINITIONS})
target_link_libraries(qcCluster PUBLIC qcMesh qcCommon gtest)

add_subdirectory(tests)