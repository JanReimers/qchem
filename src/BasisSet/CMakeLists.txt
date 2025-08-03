add_library(qcBasisSet) 
target_sources(qcBasisSet PRIVATE 
  Imp/BasisSet.C
  Internal/Imp/HeapDB.C
  Internal/Imp/Cache4.C
  Internal/Imp/ERI4.C
  Internal/Imp/BS_Common.C
  Internal/Imp/IBS_Common.C
  Internal/Imp/TIBS_Common.C
  Internal/Imp/Fit_IBS_Common.C
  Internal/Imp/DFT_IBS_Common.C
  Internal/Imp/HF_IBS_Common.C
  Internal/Imp/DHF_IBS_Common.C 
)
target_sources(qcBasisSet PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
  BasisSet.C
  Irrep_BS.C
  BasisFunction.C
  Fit_IBS.C
  HF_IBS.C
  DFT_IBS.C
  DHF_IBS.C
  Internal/Integrals.C
  Internal/ERI4.C
  Internal/Cache4.C
  Internal/IEClient.C
  Internal/IntegralEnums.C
  Internal/BS_Common.C
  Internal/IBS_Common.C
  Internal/HeapDB.C

)

set_property(TARGET qcBasisSet PROPERTY CXX_STANDARD 20)
target_link_libraries(qcBasisSet PUBLIC qcLASolver  qcSymmetry qcCluster qcMesh qcCommon) 

add_subdirectory(Atom)
add_subdirectory(Molecule)
#add_subdirectory(Lattice)

