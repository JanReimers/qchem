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
    Atom.C 
    Cluster.C 
    Molecule.C 
    UnitCell.C 
    Lattice.C 
    Internal/UnitCellMesh.C 
    Internal/AtomMesh.C 
    Internal/MoleculeMesh.C
    )
target_include_directories(qcCluster PRIVATE 
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
)
target_compile_options(qcCluster PRIVATE -Wall -Wno-unknown-pragmas 
    $<$<CONFIG:Debug>:-g -O0 >
    $<$<CONFIG:Release>: -O2 >
)
set_property(TARGET qcCluster PROPERTY CXX_STANDARD 20)
target_link_libraries(qcCluster PUBLIC qcMesh qcCommon)
