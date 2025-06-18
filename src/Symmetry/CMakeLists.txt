add_library(qcSymmetry) 
target_sources(qcSymmetry PRIVATE Atom_EC.C  ElectronCounts.C  Irrep_QNs.C  Molecule_EC.C  Okmj.C  Orbital_QNs.C  Symmetry.C  Yl.C  Ylm.C)
target_include_directories(qcSymmetry PRIVATE ${PROJECT_SOURCE_DIR}/include  ${PROJECT_SOURCE_DIR}/submodules/oml/include)
target_link_libraries(qcSymmetry PRIVATE qcCommon oml)
