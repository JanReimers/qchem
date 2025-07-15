add_library(qcCluster) 
target_sources(qcCluster PRIVATE Atom.C AtomMesh.C Molecule.C MoleculeMesh.C )
target_sources(qcCluster PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES UnitCell.C Lattice.C UnitCellMesh.C)
target_include_directories(qcCluster PRIVATE 
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/include/Imp/Cluster  
    ${PROJECT_SOURCE_DIR}/submodules/oml/include
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
)
target_compile_options(qcCluster PRIVATE -Wall -Wno-unknown-pragmas 
    $<$<CONFIG:Debug>:-g -O0>
    $<$<CONFIG:Release>: -O2>
)
set_property(TARGET qcCluster PROPERTY CXX_STANDARD 20)
target_link_libraries(qcCluster PUBLIC qcMesh qcCommon)
