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
    , itsSpin(Spin::None)
{};

SlaterExchange::SlaterExchange(double theAlpha)
    : itsAlpha(theAlpha)
    , itsSpin(Spin::None)
{};

SlaterExchange::SlaterExchange(double theAlpha, const Spin& S)
    : itsAlpha(theAlpha)
    , itsSpin(S)
{
    assert(itsSpin!=Spin::None);
};

double SlaterExchange::GetVxc(double ro) const
{
    if (itsSpin==Spin::None) ro*=0.5;
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
