enable_language(Fortran)
set(CMAKE_Fortran_FLAGS "-O2 ")
target_sources(qcAtom_BS1 PRIVATE 
  gauleg.f
  Imp/GLQuadrature.C  
  Imp/Rk.C
  Imp/SplineGrouper.C
  Imp/BS_Evaluator.C
  Imp/BS_r_Evaluator.C
  Imp/IBS_Evaluator.C
  Imp/IBS_r_Evaluator.C
)

target_sources(qcAtom_BS1 PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
  GLQuadrature.C  
  Rk.C
  SplineGrouper.C
  BS_Evaluator.C
  IBS_Evaluator.C
)

