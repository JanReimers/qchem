add_library(qcCluster) 
target_sources(qcCluster PRIVATE Atom.C AtomMesh.C Lattice.C Molecule.C MoleculeMesh.C UnitCell.C UnitCellMesh.C)
target_include_directories(qcCluster PRIVATE 
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/include/Imp/Cluster  
    ${PROJECT_SOURCE_DIR}/submodules/oml/include
    )
