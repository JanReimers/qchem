// File: LDA_XC_UT.C  Pointwise validation of hand-coded LDA exchange-correlation vs libxc.
//
// DFT-upgrade oracle groundwork (see doc + project_dft_upgrade_plan): hand-code the LDA functionals so
// we have full control, and pin them against libxc to machine precision at fixed densities.  If the
// hand-coded functionals match libxc, the XC FUNCTIONAL is correct -- so any disagreement with the NIST
// atomic LDA energies (and the historical exchange "fudge factor" in Libxc_LDA_Exchange) is NUMERICAL
// (radial grid / Vxc fit), not functional.
//
//   Exchange (Dirac/Slater, spin-unpolarized):  eps_x = -3/4 (3/pi)^{1/3} rho^{1/3},  v_x = -(3/pi)^{1/3} rho^{1/3}
//   Correlation (VWN5 paramagnetic):            Vosko-Wilk-Nusair, Can.J.Phys 58, 1200 (1980), functional V.
#include <cmath>
#include <src/xc.h>
#include "gtest/gtest.h"

namespace
{
constexpr double kPi = 3.14159265358979323846;

// --- Dirac exchange (unpolarized) ---
double DiracEps_x(double rho) { return -0.75*std::cbrt(3.0/kPi)*std::cbrt(rho); }   // energy / particle
double DiracV_x  (double rho) { return     -std::cbrt(3.0/kPi)*std::cbrt(rho); }    // potential = 4/3 eps

// --- VWN5 correlation (paramagnetic), eps_c(rho) and v_c(rho) ---
// A in Hartree (= A_Rydberg/2); b,c,x0 dimensionless.  rs = (3/4 pi rho)^{1/3}, x = sqrt(rs).
struct VWN
{
    static constexpr double A  = 0.0310907;
    static constexpr double b  = 3.72744;
    static constexpr double c  = 12.9352;
    static constexpr double x0 = -0.10498;
    static double X(double x) { return x*x + b*x + c; }

    static double Eps_c(double rho)
    {
        double rs=std::cbrt(3.0/(4.0*kPi*rho)), x=std::sqrt(rs);
        double Q=std::sqrt(4.0*c-b*b), Xx=X(x), X0=X(x0);
        double at=std::atan(Q/(2.0*x+b));
        double beta=b*x0/X0;
        double t1=std::log(x*x/Xx) + (2.0*b/Q)*at;
        double t2=std::log((x-x0)*(x-x0)/Xx) + (2.0*(b+2.0*x0)/Q)*at;
        return A*(t1 - beta*t2);
    }

    // v_c = eps_c - (x/6) d eps_c / dx   (since rs = x^2, v_c = eps_c - (rs/3) d eps_c/d rs)
    static double V_c(double rho)
    {
        double rs=std::cbrt(3.0/(4.0*kPi*rho)), x=std::sqrt(rs);
        double Q=std::sqrt(4.0*c-b*b), Xx=X(x);
        double beta=b*x0/X0_();
        // d eps_c/dx = A[ 2/x - (2x+2b)/X - beta( 2/(x-x0) - (2x+2b+2 x0)/X ) ]   (arctan terms cancel to -b/X form)
        double depsdx=A*( 2.0/x - (2.0*x+2.0*b)/Xx
                         - beta*( 2.0/(x-x0) - (2.0*x+2.0*b+2.0*x0)/Xx ) );
        return Eps_c(rho) - (x/6.0)*depsdx;
    }
    static double X0_() { return X(x0); }
};

// libxc exc & vxc for a functional id at one density (unpolarized).
void Libxc(int id, double rho, double& exc, double& vxc)
{
    xc_func_type f;
    ASSERT_EQ(xc_func_init(&f, id, XC_UNPOLARIZED), 0);
    double r[1]={rho}, e[1]={0}, v[1]={0};
    xc_lda_exc_vxc(&f, 1, r, e, v);
    exc=e[0]; vxc=v[0];
    xc_func_end(&f);
}
} // namespace

class LDA_XC : public ::testing::Test {};

// Dirac exchange must equal libxc LDA_X (id 1) to machine precision -- exchange is parameter-free.
TEST_F(LDA_XC, DiracExchangeMatchesLibxc)
{
    for (double rho : {0.01, 0.1, 0.5, 1.0, 5.0, 50.0})
    {
        double exc,vxc; Libxc(1, rho, exc, vxc);
        EXPECT_NEAR(DiracEps_x(rho), exc, 1e-12) << "rho="<<rho;
        EXPECT_NEAR(DiracV_x  (rho), vxc, 1e-12) << "rho="<<rho;
    }
}

// Hand-coded VWN5 must equal libxc LDA_C_VWN (id 7 = VWN5) to machine precision.
TEST_F(LDA_XC, VWN5CorrelationMatchesLibxc)
{
    for (double rho : {0.01, 0.1, 0.5, 1.0, 5.0, 50.0})
    {
        double exc,vxc; Libxc(7, rho, exc, vxc);
        EXPECT_NEAR(VWN::Eps_c(rho), exc, 1e-9) << "rho="<<rho;
        EXPECT_NEAR(VWN::V_c  (rho), vxc, 1e-9) << "rho="<<rho;
    }
}
