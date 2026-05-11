target_sources(qcAtom_BS1 PRIVATE 
  Imp/ExponentScaler.C
  Imp/Rk.C
  Imp/BS_Evaluator.C
  Imp/IBS_Evaluator.C
  Imp/IBS_RKBS_Evaluator.C
)

target_sources(qcAtom_BS1 PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
  ExponentScaler.C
  Integrals.C
  Rk.C
  BS_Evaluator.C
  IBS_Evaluator.C
  
)
