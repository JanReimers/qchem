// File: VWN_Correlation.C  VWN5 LDA correlation: the energy per volume, and Enzyme's derivatives of it.
//
// The functional is written once, as the value e_c(rho_up, rho_down) = rho eps_c(r_s, zeta).  The potential and the
// kernel are Enzyme's derivatives of it (the OpenWork "Enzyme REQUIRED" row, family 2, 2026-11-09): v^sigma by reverse
// mode, f^{sigma sigma'} by forward mode over that.  This unit carries the plugin (qchem_enzyme_sources in
// ../../CMakeLists.txt).  While EnzymeAD/Enzyme#3747 (a silently wrong reverse-mode gradient) is open the derivatives keep
// their finite-difference gates (tests/LDA_XC.C: v against the FD of e, f against the FD of v, f^ud = f^du, the zeta = 0
// collapse, libxc pointwise).  Its trigger, a loop-local array filled through an inner loop with a data-dependent trip
// count, is absent: the value is straight-line scalar code.
module;
#include <cmath>
#include <enzyme/enzyme>
module qchem.Hamiltonian.Internal.VWN_Correlation;

namespace qchem::Hamiltonian
{

// ---- the value ---------------------------------------------------------------------------------------------------------
// VWN5 parameters (Hartree): paramagnetic, ferromagnetic, and spin stiffness alpha_c.  G(x;A,b,c,x0) is the closed form
// shared by the three branches, with rs=(3/4 pi rho)^{1/3}, x=sqrt(rs), X(x)=x^2+b x+c.
struct VWN_Params { double A,b,c,x0; };
constexpr double kVWN_Pi = 3.14159265358979323846;
constexpr VWN_Params kVWN_P    {0.0310907 , 3.72744, 12.9352 , -0.10498   };                 // paramagnetic eps_c^P
constexpr VWN_Params kVWN_F    {0.01554535, 7.06042, 18.0578 , -0.32500   };                 // ferromagnetic eps_c^F
constexpr VWN_Params kVWN_Alpha{-1.0/(6.0*kVWN_Pi*kVWN_Pi), 1.13107, 13.0045, -0.0047584};  // spin stiffness, A=-1/(6 pi^2)
constexpr double kVWN_Den = 0.5198420997897464; // 2^{4/3}-2
constexpr double kVWN_fpp0= 1.7099209341613657; // f''(0) = 4/(9(2^{1/3}-1))

inline double VWN_X(double x, const VWN_Params& p) { return x*x + p.b*x + p.c; }
inline double VWN_G(double x, const VWN_Params& p)
{
    double Q=std::sqrt(4.0*p.c-p.b*p.b), Xx=VWN_X(x,p), X0=VWN_X(p.x0,p);
    double at=std::atan(Q/(2.0*x+p.b)), beta=p.b*p.x0/X0;
    double t1=std::log(x*x/Xx)               + (2.0*p.b/Q)*at;
    double t2=std::log((x-p.x0)*(x-p.x0)/Xx) + (2.0*(p.b+2.0*p.x0)/Q)*at;
    return p.A*(t1 - beta*t2);
}
//! x = sqrt(rs), rs = (3/4 pi rho)^{1/3}, with the constant outside the cube root: cbrt(3/(4 pi rho)) overflows to inf for
//! rho below ~1e-308 (an atom's far tail reaches the subnormals), and inf - inf in G and its derivative is NaN.
const double kVWN_rs1=std::cbrt(3.0/FourPi);   // rs at rho = 1
inline double VWN_x(double rho) { return std::sqrt(kVWN_rs1/std::cbrt(rho)); }
//! f(zeta), the spin interpolation from paramagnetic (0) to ferromagnetic (1).
inline double VWN_f(double z) { return (std::pow(1.0+z,4.0/3.0)+std::pow(1.0-z,4.0/3.0)-2.0)/kVWN_Den; }

//! eps_c^P(rho): the scalar (zeta = 0) face.
inline double VWN_EpsP(double rho) { return VWN_G(VWN_x(rho),kVWN_P); }
//! eps_c(rho, zeta), VWN Eq. 4.4 (the form libxc uses):
//! eps_c^P + alpha_c f(zeta)/f''(0) (1-zeta^4) + (eps_c^F - eps_c^P) f(zeta) zeta^4.
inline double VWN_Eps(double rho, double z)
{
    const double x=VWN_x(rho), epsP=VWN_G(x,kVWN_P), epsF=VWN_G(x,kVWN_F), ac=VWN_G(x,kVWN_Alpha);
    const double f=VWN_f(z), z4=z*z*z*z;
    return epsP + ac*(f/kVWN_fpp0)*(1.0-z4) + (epsF-epsP)*f*z4;
}

// The differentiated functions: energies PER VOLUME of physical channel densities (rho_sigma >= 0, rho > 0).  Module
// linkage, not internal: Enzyme's entry is a template declared and never defined (as in qcMesh's BeckePartition.C).
double VWN_EP(double rho)           { return rho*VWN_EpsP(rho); }
double VWN_E (double up, double dn) { const double rho=up+dn; return rho*VWN_Eps(rho,(up-dn)/rho); }

// ---- the physical point ------------------------------------------------------------------------------------------------
//! zeta = (rho_up-rho_down)/rho is a polarization FRACTION, so |zeta| > 1 is never physics: it means one channel came in
//! slightly NEGATIVE, which a band-limited rho-tilde (or any quadrature-represented density) does in the tails by ringing.
//! Fed through, it is fatal, not merely inaccurate: f(zeta) evaluates (1-zeta)^{4/3} by std::pow, NaN for a negative
//! base, and one NaN mesh point makes the whole v_c matrix non-Hermitian ("Invalid assignment to diagonal matrix element"
//! out of Blaze).  So a negative channel is mapped to the fully polarized point of the SAME total, (rho, 0) or (0, rho),
//! BEFORE differentiation: the derivative is the one at that physical point, never a derivative of the clamp (a clamp
//! inside the differentiated value would zero d zeta/d rho_sigma and drop the (1 -/+ zeta) d eps/d zeta part of v).
//! Unreachable before 2026-08-07: every polarized periodic run was fed rho_up = rho_down = rho/2 by the spin-blind
//! rho-tilde mixers (the MnO AFM collapse; doc/SymmetryUpgradePlan.md §7 step 7).
struct VWN_Point { double up, dn; };
inline VWN_Point VWN_Physical(double up, double dn)
{
    const double rho=up+dn;
    if (dn<0.0) return {rho,0.0};
    if (up<0.0) return {0.0,rho};
    return {up,dn};
}

//! The kernel's zeta threshold.  f''(zeta) ~ (1 -/+ zeta)^{-2/3} diverges at full polarization, and at EXACTLY zeta = +-1 the
//! forward-over-reverse derivative reads 0 * inf = NaN (d zeta/d rho_majority is 0 there).  The kernel therefore floors the
//! minority channel at kZetaThreshold * rho, keeping 1 -/+ zeta >= 2e-10, far above roundoff, so zeta never rounds to +-1
//! and the kernel is large but finite (libxc's zeta_threshold, with a margin).  The potential and energy need no floor.
constexpr double kZetaThreshold = 1e-10;
inline VWN_Point VWN_KernelPoint(double up, double dn)
{
    VWN_Point p=VWN_Physical(up,dn);
    const double rho=up+dn, floor=kZetaThreshold*rho;
    if (p.dn<floor) return {rho-floor,floor};
    if (p.up<floor) return {floor,rho-floor};
    return p;
}

// ---- the derivatives (Enzyme) ------------------------------------------------------------------------------------------
double VWN_VP(double rho) { return __enzyme_autodiff<double>(reinterpret_cast<void*>(VWN_EP), rho); }
double VWN_V (double up, double dn, int s)
{
    auto g=__enzyme_autodiff<enzyme::tuple<double,double>>(reinterpret_cast<void*>(VWN_E), enzyme_out, up, enzyme_out, dn);
    return s==0 ? enzyme::get<0>(g) : enzyme::get<1>(g);
}
double VWN_Vup(double up, double dn) { return VWN_V(up,dn,0); }
double VWN_Vdn(double up, double dn) { return VWN_V(up,dn,1); }

// ---- the faces -----------------------------------------------------------------------------------------------------------
double VWN_Correlation::GetVxc  (double rho) const { return rho>0.0 ? VWN_VP(rho)   : 0.0; }
double VWN_Correlation::GetEpsXc(double rho) const { return rho>0.0 ? VWN_EpsP(rho) : 0.0; }

double VWN_Correlation::GetEpsXc(double rup, double rdn, const Spin&) const
{
    if (rup==rdn) return GetEpsXc(rup+rdn);               // exact zeta = 0: the scalar path
    const double rho=rup+rdn;
    if (rho<=0.0) return 0.0;
    const VWN_Point p=VWN_Physical(rup,rdn);
    return VWN_Eps(rho,(p.up-p.dn)/rho);
}

double VWN_Correlation::GetVxc(double rup, double rdn, const Spin& s) const
{
    if (rup==rdn) return GetVxc(rup+rdn);                 // exact zeta = 0: the scalar path
    if (rup+rdn<=0.0) return 0.0;
    const VWN_Point p=VWN_Physical(rup,rdn);
    return VWN_V(p.up,p.dn,s==Spin::Down ? 1 : 0);
}

double VWN_Correlation::GetFxc(double rup, double rdn, const Spin& s, const Spin& t) const
{
    if (!(rup+rdn>1e-12)) return 0.0;                     // vanishing density: the same guard as the base class's kernel
    const VWN_Point p=VWN_KernelPoint(rup,rdn);
    void* v=reinterpret_cast<void*>(s==Spin::Down ? VWN_Vdn : VWN_Vup);
    return t==Spin::Down ? __enzyme_fwddiff<double>(v, enzyme_const, p.up, enzyme_dup, p.dn, 1.0)
                         : __enzyme_fwddiff<double>(v, enzyme_dup, p.up, 1.0, enzyme_const, p.dn);
}

} //namespace qchem::Hamiltonian
