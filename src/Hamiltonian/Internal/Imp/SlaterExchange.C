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

std::ostream& SlaterExchange::Write(std::ostream& os) const
{
    os << itsAlpha << " ";
    return os;
}


} //namespace
