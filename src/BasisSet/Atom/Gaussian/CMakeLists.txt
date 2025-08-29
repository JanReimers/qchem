target_sources(qcAtom_BS PRIVATE 
  Imp/ExponentScaler.C
  Imp/Rk.C
  NR/Imp/BS_Evaluator.C
  NR/Imp/IBS_Evaluator.C
  NR/Imp/Gaussian_BS.C
  NR/Imp/Gaussian_IBS.C
  RKB/Imp/Gaussian_BS.C
  RKB/Imp/Gaussian_IBS.C
)

target_sources(qcAtom_BS PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
  ExponentScaler.C
  Integrals.C
  Rk.C
  NR/BS_Evaluator.C
  NR/IBS_Evaluator.C
  NR/Gaussian_BS.C
  RKB/Gaussian_BS.C

)