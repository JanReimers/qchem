target_sources(qcAtom_BS PRIVATE 
  Imp/ExponentScaler.C
  Imp/Rk.C
  Imp/BS_Evaluator.C
  Imp/IBS_Evaluator.C
  Imp/IBS_RKBS_Evaluator.C
)

target_sources(qcAtom_BS PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
  ExponentScaler.C
  Integrals.C
  Rk.C
  BS_Evaluator.C
  IBS_Evaluator.C
  IBS_RKBS_Evaluator.C
  
)
target_compile_options(qchem PRIVATE 
    $<$<CONFIG:Debug>:-g -O0>
    $<$<CONFIG:Release>: -O2 -pg>
)