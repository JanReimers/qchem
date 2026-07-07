add_library(qcBasisSet) 
target_sources(qcBasisSet PRIVATE 
  Imp/BasisSet.C
  Imp/IrrepBasisSet.C
  Imp/Orbital_1E_IBS.C
  Imp/Orbital_HF_IBS.C
  Imp/Orbital_DFT_IBS.C
  Imp/Orbital_DHF_IBS.C
  Imp/Fit_IBS.C
  Imp/SymmetryAdapted_IBS.C
  Internal/Imp/Orbital_DHF_IBS.C
  Internal/Imp/DB_Cache.C
  Internal/Imp/Cache4.C
  Internal/Imp/Cache2.C
  Internal/Imp/Cache3.C
  Internal/Imp/ERI3.C
  Internal/Imp/ERI4.C
  Internal/Imp/ERI4T.C
  Internal/Imp/DB_Cache_RAM.C
  Internal/Imp/DB_Cache_RAM_Report.C

)
target_sources(qcBasisSet PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
  BasisSet.C
  IrrepBasisSet.C
  Orbital_1E_IBS.C
  Band_DFT_IBS.C
  Band_FT_IBS.C
  G_FieldEvaluator.C
  Orbital_HF_IBS.C
  Orbital_DFT_IBS.C
  Orbital_DHF_IBS.C
  Fit_IBS.C
  SymmetryAdapted_IBS.C
  Internal/BasisSetImp.C
  Internal/IrrepBasisSetImp.C
  Internal/Orbital_DHF_IBS.C
  Internal/Cache4.C
  Internal/Cache2.C
  Internal/Cache3.C
  Internal/ERI3.C
  Internal/GMap.C
  Internal/ERI4.C
  Internal/ERI4T.C
  Internal/IntegralEnums.C
  Internal/DB_Cache.C
  Internal/DB_Cache_RAM.C
  )
target_include_directories(qcBasisSet PRIVATE ${PROJECT_SOURCE_DIR}/src)
set_property(TARGET qcBasisSet PROPERTY CXX_STANDARD 20)
target_include_directories(qcBasisSet PRIVATE ${PROJECT_SOURCE_DIR}/submodules/blaze)
target_link_libraries(qcBasisSet PUBLIC qcElConfig qcStructure qcMesh)
target_compile_options(qcBasisSet PRIVATE ${CXX_COMMON_COMPILE_OPTIONS})
target_compile_definitions(qcBasisSet PUBLIC ${CXX_COMMON_COMPILE_DEFINITIONS})
add_subdirectory(Atom)
add_subdirectory(Molecule)
add_subdirectory(Lattice_3D)

