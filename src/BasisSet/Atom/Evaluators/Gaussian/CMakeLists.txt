target_sources(qcAtom_BS PRIVATE 
  Internal/Imp/ExponentScaler.C
  Internal/Imp/Rk.C
  Imp/Evaluator.C
  Imp/RKBS_Evaluator.C
)

target_sources(qcAtom_BS PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
  Internal/ExponentScaler.C
  Internal/Integrals.C
  Internal/Rk.C
  Evaluator.C
  
)
