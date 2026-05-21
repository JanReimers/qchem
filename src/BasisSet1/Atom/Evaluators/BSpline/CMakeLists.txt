enable_language(Fortran)
set(CMAKE_Fortran_FLAGS "-O2 ")
target_sources(qcAtom_BS1 PRIVATE 
  Internal/gauleg.f
  Internal/Imp/GLQuadrature1.C  
  Internal/Imp/Rk.C
  Internal/Imp/SplineGrouper.C
  Imp/IBS_Evaluator.C
  Imp/IBS_r_Evaluator.C
)

target_sources(qcAtom_BS1 PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
  Internal/GLQuadrature.C  
  Internal/Rk.C
  Internal/SplineGrouper.C
  IBS_Evaluator.C
)

