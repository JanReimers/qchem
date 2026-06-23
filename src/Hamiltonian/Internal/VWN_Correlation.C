// File: VWN_Correlation.C  Vosko-Wilk-Nusair (VWN5) LDA correlation functional.
module;
#include <cmath>
#include <ostream>
export module qchem.Hamiltonian.Internal.VWN_Correlation;
import qchem.Hamiltonian.Internal.ExFunctional;
import qchem.ChargeDensity;
import qchem.Vector3D;
import qchem.Math;   // Pi, FourPi

export namespace qchem::Hamiltonian
{

//! \brief VWN5 (Vosko-Wilk-Nusair, Can.J.Phys 58, 1200 (1980), functional V) LDA CORRELATION, paramagnetic
//! (spin-unpolarized).  Validated pointwise against libxc LDA_C_VWN (see LDA_XC_UT).  Provides BOTH the
//! potential \f$v_c\f$ (GetVxc) and the energy density \f$\varepsilon_c\f$ (GetEpsXc) -- the latter is
//! essential because \f$E_c=\int\varepsilon_c\rho\f$ and \f$\varepsilon_c\neq\tfrac34 v_c\f$, so the
//! exchange-virial energy shortcut used for Slater exchange is invalid for correlation.
class VWN_Correlation : public ExFunctional
{
    typedef Vector3D<double> Vec3;
public:
    VWN_Correlation() {}

    virtual double operator()(const Vec3& r) const { return GetVxc((*itsChargeDensity)(r)); }
    virtual Vec3   Gradient  (const Vec3&  ) const { return Vec3(0,0,0); } // unused: the fit uses values

    virtual double GetVxc  (double rho) const { return rho>0.0 ? Vc (rho) : 0.0; }
    virtual double GetEpsXc(double rho) const { return rho>0.0 ? Eps(rho) : 0.0; }

    virtual std::ostream& Write(std::ostream& os) const { return os << "VWN5"; }

private:
    // VWN5 paramagnetic parameters (Hartree).  rs=(3/4 pi rho)^{1/3}, x=sqrt(rs), X(x)=x^2+b x+c.
    static constexpr double A=0.0310907, b=3.72744, c=12.9352, x0=-0.10498;
    static double X(double x) { return x*x + b*x + c; }

    static double Eps(double rho)   // eps_c
    {
        double rs=std::cbrt(3.0/(FourPi*rho)), x=std::sqrt(rs);
        double Q=std::sqrt(4.0*c-b*b), Xx=X(x), at=std::atan(Q/(2.0*x+b)), beta=b*x0/X(x0);
        double t1=std::log(x*x/Xx)            + (2.0*b/Q)*at;
        double t2=std::log((x-x0)*(x-x0)/Xx)  + (2.0*(b+2.0*x0)/Q)*at;
        return A*(t1 - beta*t2);
    }
    static double Vc(double rho)    // v_c = eps_c - (x/6) d eps_c/dx   (rs = x^2)
    {
        double rs=std::cbrt(3.0/(FourPi*rho)), x=std::sqrt(rs), Xx=X(x), beta=b*x0/X(x0);
        double depsdx=A*( 2.0/x - (2.0*x+2.0*b)/Xx - beta*(2.0/(x-x0) - (2.0*x+2.0*b+2.0*x0)/Xx) );
        return Eps(rho) - (x/6.0)*depsdx;
    }
};

} //namespace
