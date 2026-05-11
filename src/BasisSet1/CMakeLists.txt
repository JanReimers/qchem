add_library(qcBasisSet1) 
target_sources(qcBasisSet1 PRIVATE 
  Imp/DB_Cache.C
  Imp/BasisSet.C
  Imp/IrrepBasisSet.C
  Imp/Orbital_1E_IBS.C
  Imp/Orbital_HF_IBS.C
  Imp/Orbital_DFT_IBS.C
  Imp/Orbital_DHF_IBS.C
  Imp/Fit_IBS.C
  Internal/Imp/Orbital_DHF_IBS.C
  Internal/Imp/Cache4.C
  Internal/Imp/ERI3.C
  Internal/Imp/ERI4.C
  Internal/Imp/ERI4T.C

)
target_sources(qcBasisSet1 PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
  DB_Cache.C
  BasisSet.C
  IrrepBasisSet.C
  Orbital_1E_IBS.C
  Orbital_HF_IBS.C
  Orbital_DFT_IBS.C
  Orbital_DHF_IBS.C
  Fit_IBS.C
  Internal/BasisSetImp.C
  Internal/IrrepBasisSetImp.C
  Internal/Orbital_DHF_IBS.C
  Internal/Cache4.C
  Internal/ERI3.C
  Internal/ERI4.C
  Internal/ERI4T.C
  Internal/IntegralEnums.C
  )
target_include_directories(qcBasisSet1 PRIVATE ${PROJECT_SOURCE_DIR}/src)
set_property(TARGET qcBasisSet1 PROPERTY CXX_STANDARD 20)
target_include_directories(qcBasisSet1 PRIVATE ${PROJECT_SOURCE_DIR}/submodules/blaze)
target_link_libraries(qcBasisSet1 PUBLIC qcLASolver qcSymmetry qcCluster qcMesh qcCommon) 
target_compile_options(qchem PRIVATE 
    $<$<CONFIG:Debug>:-g -O0>
    $<$<CONFIG:Release>: -O2 -pg>
)
add_subdirectory(Atom)
add_subdirectory(Molecule)
# add_subdirectory(Lattice)

