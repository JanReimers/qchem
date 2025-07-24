target_sources(qcAtom_BS PRIVATE 
Imp/Slater_BS.C
Imp/Slater_IBS.C
Imp/Slater_BF.C 
Imp/Gaussian_IBS.C 
Imp/Gaussian_BS.C
Imp/Gaussian_BF.C 
Imp/BSpline_BS.C
Imp/BSpline_IBS.C
Imp/BSpline_BF.C 
)

target_sources(qcAtom_BS PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
  IE_HF_Angular.C
  Gaussian_BS.C
  Slater_BS.C
  BSpline_BS.C
  
)