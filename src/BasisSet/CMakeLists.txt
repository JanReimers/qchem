add_library(qcBasisSet) 
target_sources(qcBasisSet PRIVATE 
  Imp/BasisSet.C
  Imp/IrrepBasisSet.C
  Imp/DB_Cache1.C
  Imp/Orbital_HF_IBS1.C
  Internal/Imp/HeapDB.C
  Internal/Imp/DB_Cache.C
  Internal/Imp/Cache4.C
  Internal/Imp/ERI4.C
  Internal/Imp/ERI4T.C
  Internal/Imp/BasisSet.C
  Internal/Imp/BasisSet1.C
  Internal/Imp/IrrepBasisSet.C
  Internal/Imp/Fit_IBS_Common.C
  Internal/Imp/DFT_IBS_Common.C
  Internal/Imp/HF_IBS_Common.C
  Internal/Imp/DHF_IBS_Common.C
)
target_sources(qcBasisSet PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
  BasisSet.C
  IrrepBasisSet.C
  Fit_IBS.C
  Orbital_HF_IBS.C
  Orbital_DFT_IBS.C
  Orbital_DHF_IBS.C
  Orbital_1E_IBS.C
  DB_Cache1.C
  BasisSet1.C
  IrrepBasisSet1.C
  Orbital_1E_IBS1.C
  Orbital_HF_IBS1.C
  Internal/ERI3.C
  Internal/ERI4T.C
  Internal/ERI4.C
  Internal/Cache4.C
  Internal/IntegralEnums.C
  Internal/BS_Common.C
  Internal/BS_Common1.C
  Internal/IrrepBasisSet.C
  Internal/IrrepBasisSet1.C
  Internal/DB_Cache.C
  Internal/HeapDB.C
)
target_include_directories(qcBasisSet PRIVATE ${PROJECT_SOURCE_DIR}/src)
set_property(TARGET qcBasisSet PROPERTY CXX_STANDARD 20)
target_include_directories(qcBasisSet PRIVATE ${PROJECT_SOURCE_DIR}/submodules/blaze)
target_link_libraries(qcBasisSet PUBLIC qcLASolver qcSymmetry qcCluster qcMesh qcCommon) 
target_compile_options(qchem PRIVATE 
    $<$<CONFIG:Debug>:-g -O0>
    $<$<CONFIG:Release>: -O2 -pg>
)
add_subdirectory(Atom)
add_subdirectory(Molecule)
add_subdirectory(Lattice)

