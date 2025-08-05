target_sources(qcAtom_BS PRIVATE 
    # Imp/Gaussian_BS.C
    # Imp/Gaussian_IBS.C
    #Imp/BSpline_BS.C
    #Imp/BSpline_IBS.C
)
target_sources(qcAtom_BS PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
    IE_HF_Angular.C
    # Gaussian_BS.C
    #BSpline_BS.C
)