enable_language(Fortran)
set(CMAKE_Fortran_FLAGS "-O2 ")
target_sources(qcAtom_BS PRIVATE 
  Internal/gauleg.f
  Internal/Imp/GLQuadrature1.C  
  Internal/Imp/Rk.C
  Internal/Imp/SplineGrouper.C
  Internal/Imp/EvaluatorCommon.C
  Imp/Evaluator.C
  Imp/Evaluator_r.C
)

target_sources(qcAtom_BS PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
  Internal/GLQuadrature.C  
  Internal/Rk.C
  Internal/SplineGrouper.C
  Internal/EvaluatorCommon.C
  Evaluator.C
  Evaluator_r.C
)

