target_sources(qcAtom_BS PRIVATE 
    Imp/Slater_IBS.C
    Imp/Gaussian_BS.C
    Imp/Slater_BS.C
    Imp/Slater_IE.C
    Imp/Gaussian_IE.C
    Imp/Gaussian_IBS.C
)
target_sources(qcAtom_BS PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
    Gaussian_BS.C
    Slater_BS.C
)