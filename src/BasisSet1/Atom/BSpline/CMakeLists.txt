enable_language(Fortran)
set(CMAKE_Fortran_FLAGS "-O2 -askjdfhalsdjf")
target_sources(qcAtom_BS1 PRIVATE 
 
)

target_sources(qcAtom_BS1 PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
  
  NR/BSpline_BS1.C
)

