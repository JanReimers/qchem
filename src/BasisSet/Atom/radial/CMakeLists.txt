target_sources(qcAtom_BS PRIVATE 
    Slater/Imp/Rk.C
    Slater/Imp/ExponentScaler.C
    Gaussian/Imp/Rk.C
    Gaussian/Imp/ExponentScaler.C

    BSpline/GLQuadrature.C
    BSpline/gauleg.f
    BSpline/Imp/Rk.C
    BSpline/Imp/SplineGrouper.C
    BSpline/Imp/GLQuadrature.C

    Imp/Slater_IBS.C
    Imp/Slater_BS.C
    Imp/Gaussian_IBS.C
    Imp/Gaussian_BS.C
    Imp/BSpline_IBS.C
    Imp/BSpline_BS.C
    Imp/ExponentGrouper.C
    Imp/Exponential_Evaluator.C
)
target_sources(qcAtom_BS PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
    PascalTriangle.C
    FillPower.C
    Gaussian/Rk.C
    Gaussian/ExponentScaler.C
    Gaussian/Integrals.C
    Slater/Rk.C
    Slater/ExponentScaler.C
    Slater/Integrals.C
    BSpline/Rk.C
    BSpline/SplineGrouper.C
    BSpline/GLQuadrature.C

    Slater_IBS.C
    Slater_BS.C
    Gaussian_IBS.C
    Gaussian_BS.C
    BSpline_IBS.C
    BSpline_BS.C
    ExponentGrouper.C
    Exponential_Evaluator.C
)