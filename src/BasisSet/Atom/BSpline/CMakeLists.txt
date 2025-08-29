target_sources(qcAtom_BS PRIVATE 
  gauleg.f
  Imp/GLQuadrature.C  
  Imp/Rk.C
  Imp/SplineGrouper.C
  NR/Imp/BS_Evaluator.C
  NR/Imp/IBS_Evaluator.C
  NR/Imp/BSpline_BS.C
  NR/Imp/BSpline_IBS.C
)

target_sources(qcAtom_BS PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
  GLQuadrature.C  
  Rk.C
  SplineGrouper.C
  NR/BS_Evaluator.C
  NR/IBS_Evaluator.C
  NR/BSpline_BS.C
)