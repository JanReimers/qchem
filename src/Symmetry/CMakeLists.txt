add_library(qc_Symmetry) 
target_sources(qc_Symmetry PRIVATE Atom_EC.C  ElectronCounts.C  Irrep_QNs.C  Molecule_EC.C  Okmj.C  Orbital_QNs.C  Symmetry.C  Yl.C  Ylm.C)
target_include_directories(qc_Symmetry PRIVATE ${PROJECT_SOURCE_DIR}/include  ${PROJECT_SOURCE_DIR}/submodules/oml/include)

