target_sources(qcAtom_BS PRIVATE 
  Imp/ExponentScaler.C
  Imp/Rk.C
  NR/Imp/BS_Evaluator.C
  NR/Imp/IBS_Evaluator.C
  NR/Imp/Slater_BS.C
  NR/Imp/Slater_IBS.C
  RKB/Imp/Slater_BS.C
  RKB/Imp/Slater_IBS.C
)

target_sources(qcAtom_BS PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
  ExponentScaler.C
  Integrals.C
  Rk.C
  NR/BS_Evaluator.C
  NR/IBS_Evaluator.C
  NR/Slater_BS.C
  RKB/Slater_BS.C

)