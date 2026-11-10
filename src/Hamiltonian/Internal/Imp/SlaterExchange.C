// File: SlaterExchange.C  Slater exchange: the energy, and Enzyme's derivatives of it.
//
// The functional is written once, as eps_x(rho); the potential d(rho eps_x)/d rho is Enzyme's reverse mode and the kernel
// forward mode over that (the OpenWork "Enzyme REQUIRED" row, family 2, 2026-11-09).  This unit carries the plugin
// (qchem_enzyme_sources in ../../CMakeLists.txt).  The FD gates stay while EnzymeAD/Enzyme#3747 is open (tests/LDA_XC.C);
// its trigger (a loop-local array filled through an inner loop) is absent from this straight-line code.
module;
#include <cmath>
#include <iostream>
#include <enzyme/enzyme>
module qchem.Hamiltonian.Internal.SlaterExchange;

namespace qchem::Hamiltonian
{

SlaterExchange::SlaterExchange()
    : itsAlpha(0)
{};

SlaterExchange::SlaterExchange(double theAlpha)
    : itsAlpha(theAlpha)
{};

//! eps_x(rho) = -(9/4) alpha (3 rho/8 pi)^{1/3}: the closed-shell energy per particle (each channel holds rho/2).
const double kSlaterC=std::cbrt(3.0/(2.0*FourPi));
inline double SlaterEps(double rho, double alpha) { return -2.25*alpha*kSlaterC*std::cbrt(rho); }
// The differentiated functions.  Module linkage, not internal: Enzyme's entry is a template declared and never defined.
//! e_x = rho eps_x = -(9/4) alpha (3/8 pi)^{1/3} rho^{4/3}, written as ONE power so neither derivative divides by rho: Enzyme takes
//! d x^y/dx = y x^{y-1}, finite down to the smallest subnormal, which an atom's far tail reaches (1e-323, the H atom on an SG mesh).
//! As rho cbrt(3 rho/8 pi) the product underflowed to 0 there and v read cbrt(x)/(3x) = 0/0, and f's 1/rho^2 overflowed to inf - inf.
double SlaterE(double rho, double alpha) { return -2.25*alpha*kSlaterC*std::pow(rho,4.0/3.0); }
double SlaterV(double rho, double alpha)
{
    return __enzyme_autodiff<double>(reinterpret_cast<void*>(SlaterE), rho, enzyme_const, alpha);
}

double SlaterExchange::GetVxc(double ro) const
{
    return ro>0.0 ? SlaterV(ro,itsAlpha) : 0.0;
}

double SlaterExchange::GetEpsXc(double ro) const
{
    return ro>0.0 ? SlaterEps(ro,itsAlpha) : 0.0;
}

double SlaterExchange::GetFxc(double up, double dn, const Spin& s, const Spin& t) const
{
    if (s!=t) return 0.0;                              // channel-separable: no cross-spin kernel
    const double rs=(s==Spin::Down ? dn : up);
    if (!(rs>0.0)) return 0.0;
    // f = d/drho_s [ v_x(2 rho_s) ] = 2 v_x'(2 rho_s)
    return 2.0*__enzyme_fwddiff<double>(reinterpret_cast<void*>(SlaterV), enzyme_dup, 2.0*rs, 1.0, enzyme_const, itsAlpha);
}

std::ostream& SlaterExchange::Write(std::ostream& os) const
{
    os << itsAlpha << " ";
    return os;
}


} //namespace
