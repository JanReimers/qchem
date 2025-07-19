add_library(qcBasisSet) 
target_sources(qcBasisSet PRIVATE 
  DataBase/HeapDBImp.C
  Imp/Cache4Imp.C
  Imp/ERI4Imp.C
  Imp/BasisSet.C
  Common/Imp/BS_Common.C
  Common/Imp/IBS_Common.C
  Common/Imp/TIBS_Common.C
  Common/Imp/Fit_IBS_Common.C
  Common/Imp/DFT_IBS_Common.C
  Common/Imp/HF_IBS_Common.C
  Common/Imp/DHF_IBS_Common.C
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
  Imp/Cache4.C
  Imp/IEClient.C
  Imp/IntegralEnums.C
  Common/BS_Common.C
  Common/IBS_Common.C
  DataBase/HeapDB.C

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

