add_library(qcBasisSet1) 
target_sources(qcBasisSet1 PRIVATE 
  Imp/DB_Cache1.C
  Imp/Orbital_HF_IBS1.C
)
target_sources(qcBasisSet1 PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
  DB_Cache1.C
  BasisSet1.C
  IrrepBasisSet1.C
  Orbital_1E_IBS1.C
  Orbital_HF_IBS1.C
  Internal/BS_Common1.C
  # Internal/IrrepBasisSet1.C
)
target_include_directories(qcBasisSet1 PRIVATE ${PROJECT_SOURCE_DIR}/src)
set_property(TARGET qcBasisSet1 PROPERTY CXX_STANDARD 20)
target_include_directories(qcBasisSet1 PRIVATE ${PROJECT_SOURCE_DIR}/submodules/blaze)
target_link_libraries(qcBasisSet1 PUBLIC qcBasisSet qcLASolver qcSymmetry qcCluster qcMesh qcCommon) 
target_compile_options(qchem PRIVATE 
    $<$<CONFIG:Debug>:-g -O0>
    $<$<CONFIG:Release>: -O2 -pg>
)
add_subdirectory(Atom)
# add_subdirectory(Molecule)
# add_subdirectory(Lattice)

