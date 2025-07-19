add_library(qcBasisSet) 
target_sources(qcBasisSet PRIVATE 
  BasisSet.C
      BS_Common.C
      IBS_Common.C
    TIBS_Common.C
  HF_IBS_Common.C
  DHF_IBS_Common.C
  DFT_IBS_Common.C
  Fit_IBS_Common.C
  DataBase/HeapDB.C
  Cache4.C
  Imp/ERI4Imp.C
  Imp/BasisSet.C

)
target_sources(qcBasisSet PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
  BasisSet.C
  Irrep_BS.C
  BasisFunction.C
  Fit_IBS.C
  HF_IBS.C
  DFT_IBS.C
  DHF_IBS.C
  Integrals.C
  Imp/ERI4.C
  Imp/IntegralEnums.C
)

target_include_directories(qcBasisSet PRIVATE 
    ${PROJECT_SOURCE_DIR}/include  
    ${PROJECT_SOURCE_DIR}/include/Imp/BasisSet 
    ${PROJECT_SOURCE_DIR}/submodules/oml/include 
    ${PROJECT_SOURCE_DIR}/submodules/BSplinebasis/include
    ${PROJECT_SOURCE_DIR}/submodules/wignerSymbols/include
    ${PROJECT_SOURCE_DIR}/submodules/json/single_include
)
set_property(TARGET qcBasisSet PROPERTY CXX_STANDARD 20)
target_link_libraries(qcBasisSet PUBLIC qcSymmetry qcCluster qcMesh qcCommon) 

add_subdirectory(Atom)
add_subdirectory(Molecule)

