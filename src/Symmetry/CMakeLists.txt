add_library(qcSymmetry) 
target_sources(qcSymmetry PRIVATE Atom_EC.C  ElectronCounts.C  Irrep_QNs.C  Molecule_EC.C  Okmj.C  Orbital_QNs.C  Symmetry.C  Yl.C  Ylm.C)
target_include_directories(qcSymmetry PRIVATE ${PROJECT_SOURCE_DIR}/include  ${PROJECT_SOURCE_DIR}/submodules/oml/include)
target_link_libraries(qcSymmetry PRIVATE qcCommon omlm)
set_property(TARGET qcSymmetry PROPERTY CXX_STANDARD 20)
# This linker association is required in order for cmake/gcc to find .gcm binaru module  interfact (BMI) files!!!!!
# target_link_libraries(qcMesh PUBLIC qcMesh qcCommon) 
