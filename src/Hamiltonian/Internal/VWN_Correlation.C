// File: VWN_Correlation.C  Vosko-Wilk-Nusair (VWN5) LDA correlation functional.
module;
#include <cmath>
#include <ostream>
export module qchem.Hamiltonian.Internal.VWN_Correlation;
import qchem.Hamiltonian.Internal.ExFunctional;
import qchem.ChargeDensity;
import qchem.Vector3D;
import qchem.Symmetry.Spin;   // Spin (the requested correlation channel)
import qchem.Math;   // Pi, FourPi

export namespace qchem::Hamiltonian
{

//! \brief VWN5 (Vosko-Wilk-Nusair, Can.J.Phys 58, 1200 (1980), functional V) LDA CORRELATION.
//!
//! SPIN-NATIVE: the primary formulation is the two-channel \f$\varepsilon_c(\rho_\uparrow,\rho_\downarrow)\f$
//! with the per-channel potentials \f$v_c^\sigma\f$ (GetEpsC / GetVc).  The spin-unpolarized scalar face
//! (GetVxc / GetEpsXc) is the \f$\zeta=0\f$ collapse and is kept BYTE-IDENTICAL to the historical
//! paramagnetic implementation (it drives the unpolarized \c Ham_DFTcorr_U SCF anchors).
//!
//! Validated pointwise against libxc LDA_C_VWN (see LDA_XC_UT), both XC_UNPOLARIZED and XC_POLARIZED.
//! Provides BOTH the potential \f$v_c\f$ and the energy density \f$\varepsilon_c\f$ -- the latter is
//! essential because \f$E_c=\int\varepsilon_c\rho\f$ and \f$\varepsilon_c\neq\tfrac34 v_c\f$, so the
//! exchange-virial energy shortcut used for Slater exchange is invalid for correlation.
//!
//! Spin interpolation (VWN Eq. 4.4, the form used by libxc):
//! \f[ \varepsilon_c(r_s,\zeta)=\varepsilon_c^P + \alpha_c(r_s)\frac{f(\zeta)}{f''(0)}(1-\zeta^4)
//!                              + [\varepsilon_c^F-\varepsilon_c^P]\,f(\zeta)\,\zeta^4 \f]
//! with \f$f(\zeta)=\frac{(1+\zeta)^{4/3}+(1-\zeta)^{4/3}-2}{2^{4/3}-2}\f$ interpolating paramagnetic
//! (\f$\zeta=0\f$) to ferromagnetic (\f$\zeta=1\f$), and \f$\alpha_c\f$ the spin stiffness.
class VWN_Correlation : public ExFunctional, public virtual SpinCorrelation
{
    typedef Vector3D<double> Vec3;
public:
    VWN_Correlation() {}

    virtual double operator()(const Vec3& r) const { return GetVxc((*itsChargeDensity)(r)); }
    virtual Vec3   Gradient  (const Vec3&  ) const { return Vec3(0,0,0); } // unused: the fit uses values

    // --- spin-unpolarized scalar face (zeta=0 collapse; BYTE-IDENTICAL to the historical code) ---
    virtual double GetVxc  (double rho) const { return rho>0.0 ? Vc (rho) : 0.0; }
    virtual double GetEpsXc(double rho) const { return rho>0.0 ? Eps(rho) : 0.0; }

    // --- spin-native two-channel face (the primary formulation; SpinCorrelation) ---
    //! Correlation energy density \f$\varepsilon_c(\rho_\uparrow,\rho_\downarrow)\f$.
    virtual double GetEpsC(double rup, double rdn) const
    {
        double rho=rup+rdn;
        return rho>0.0 ? EvalRZ(rho,Zeta(rup,rdn)).eps : 0.0;
    }
    //! Channel correlation potential \f$v_c^\sigma=\varepsilon_c+\rho\,\partial\varepsilon_c/\partial\rho_\sigma\f$.
    //! Note v_c^sigma COUPLES both channels (through r_s and zeta) -- it is NOT a function of rho_sigma alone.
    virtual double GetVc(double rup, double rdn, const Spin& s) const
    {
        double rho=rup+rdn;
        if (rho<=0.0) return 0.0;
        double z=Zeta(rup,rdn);
        Eval e=EvalRZ(rho,z);
        double rs=std::cbrt(3.0/(FourPi*rho));
        double vc=e.eps - (rs/3.0)*e.depsdrs;          // the rs-derivative (spin-common) part
        // + the zeta-derivative part:  rho d zeta/d rho_up = (1-zeta),  rho d zeta/d rho_down = -(1+zeta)
        return s==Spin::Up ? vc + (1.0-z)*e.depsdz
                           : vc - (1.0+z)*e.depsdz;
    }

    virtual std::ostream& Write(std::ostream& os) const { return os << "VWN5"; }

private:
    //! One VWN function G(x;A,b,c,x0) -- the closed form shared by the para/ferro/stiffness branches.
    struct Params { double A,b,c,x0; };

    // VWN5 parameters (Hartree): paramagnetic, ferromagnetic, and spin stiffness alpha_c.
    static constexpr double kPi = 3.14159265358979323846;
    static constexpr Params P    {0.0310907 , 3.72744, 12.9352 , -0.10498   }; // paramagnetic eps_c^P
    static constexpr Params F    {0.01554535, 7.06042, 18.0578 , -0.32500   }; // ferromagnetic eps_c^F
    static constexpr Params Alpha{-1.0/(6.0*kPi*kPi), 1.13107, 13.0045, -0.0047584}; // spin stiffness, A=-1/(6 pi^2)

    // f(zeta) spin-interpolation function and its derivative; f''(0) sets the stiffness normalization.
    static double fz  (double z) { return (std::pow(1.0+z,4.0/3.0)+std::pow(1.0-z,4.0/3.0)-2.0)/kDen; }
    static double dfz (double z) { return (4.0/3.0)*(std::cbrt(1.0+z)-std::cbrt(1.0-z))/kDen; }
    static constexpr double kDen = 0.5198420997897464; // 2^{4/3}-2
    static constexpr double kfpp0= 1.7099209341613657; // f''(0) = 4/(9(2^{1/3}-1))

    static double Zeta(double rup,double rdn) { double rho=rup+rdn; return rho>0.0 ? (rup-rdn)/rho : 0.0; }

    // G(x;p) and dG/dx for one parameter set.  rs=(3/4 pi rho)^{1/3}, x=sqrt(rs), X(x)=x^2+b x+c.
    static double Xof(double x, const Params& p) { return x*x + p.b*x + p.c; } // X(x)=x^2+b x+c (exact op order)
    static double Gval(double x, const Params& p)
    {
        double Q=std::sqrt(4.0*p.c-p.b*p.b), Xx=Xof(x,p), X0=Xof(p.x0,p);
        double at=std::atan(Q/(2.0*x+p.b)), beta=p.b*p.x0/X0;
        double t1=std::log(x*x/Xx)               + (2.0*p.b/Q)*at;
        double t2=std::log((x-p.x0)*(x-p.x0)/Xx) + (2.0*(p.b+2.0*p.x0)/Q)*at;
        return p.A*(t1 - beta*t2);
    }
    static double dGdx(double x, const Params& p)
    {
        double Xx=Xof(x,p), X0=Xof(p.x0,p), beta=p.b*p.x0/X0;
        return p.A*( 2.0/x - (2.0*x+2.0*p.b)/Xx - beta*(2.0/(x-p.x0) - (2.0*x+2.0*p.b+2.0*p.x0)/Xx) );
    }

    // --- the scalar (paramagnetic) face: G(x;P), kept byte-identical to the historical Eps/Vc ---
    static double Eps(double rho)        // eps_c^P
    {
        double x=std::sqrt(std::cbrt(3.0/(FourPi*rho)));
        return Gval(x,P);
    }
    static double Vc(double rho)         // v_c^P = eps_c^P - (x/6) d eps_c^P/dx   (rs = x^2)
    {
        double x=std::sqrt(std::cbrt(3.0/(FourPi*rho)));
        return Gval(x,P) - (x/6.0)*dGdx(x,P);
    }

    // --- the spin-native evaluation: eps_c and its (r_s, zeta) partials at one point ---
    struct Eval { double eps, depsdrs, depsdz; };
    static Eval EvalRZ(double rho, double z)
    {
        double rs=std::cbrt(3.0/(FourPi*rho)), x=std::sqrt(rs);
        double epsP=Gval(x,P), epsF=Gval(x,F), ac=Gval(x,Alpha);
        // d/d rs = (1/2x) d/dx
        double dP=dGdx(x,P)/(2.0*x), dF=dGdx(x,F)/(2.0*x), dA=dGdx(x,Alpha)/(2.0*x);
        double f=fz(z), df=dfz(z), z3=z*z*z, z4=z3*z;
        double g  = (f/kfpp0)*(1.0-z4);                       // alpha_c weight
        double dg = (df/kfpp0)*(1.0-z4) + (f/kfpp0)*(-4.0*z3);// d/d zeta of that weight
        double h  = f*z4;                                     // (epsF-epsP) weight
        double dh = df*z4 + f*4.0*z3;                         // d/d zeta of that weight
        Eval e;
        e.eps     = epsP + ac*g  + (epsF-epsP)*h;
        e.depsdrs = dP   + dA*g  + (dF  -dP )*h;
        e.depsdz  =        ac*dg + (epsF-epsP)*dh;
        return e;
    }
};

} //namespace
