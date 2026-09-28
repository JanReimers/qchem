// File: SlaterExchange.C  Slater exchange potential.
module;
#include <cassert>
#include <iostream>
module qchem.Hamiltonian.Internal.SlaterExchange;
import qchem.Math;

namespace qchem::Hamiltonian
{

SlaterExchange::SlaterExchange()
    : itsAlpha(0)
{};

SlaterExchange::SlaterExchange(double theAlpha)
    : itsAlpha(theAlpha)
{};

double SlaterExchange::GetVxc(double ro) const
{
    ro*=0.5;                    // the closed-shell face: each channel carries half the total
    double ret=0;
    if (ro > 0.0)
    {
        ret=-3.0 * itsAlpha * pow(3.0*ro/FourPi , 1.0/3.0);
    }
    return ret;
}

double SlaterExchange::GetFxc(double up, double dn, const Spin& s, const Spin& t) const
{
    if (s!=t) return 0.0;                              // channel-separable: no cross-spin kernel
    const double rs=(s==Spin::Down ? dn : up);
    if (!(rs>0.0)) return 0.0;
    return GetVxc(2.0*rs)/(3.0*rs);                    // d/drho_s [ v_x(2 rho_s) ],  v_x ~ rho^{1/3}
}

std::ostream& SlaterExchange::Write(std::ostream& os) const
{
    os << itsAlpha << " ";
    return os;
}


} //namespace
