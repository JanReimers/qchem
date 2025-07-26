// File: Imp/ExchangeFunctional.C   Exchange potential for DFT.
module;
#include <cassert>
module qchem.Hamiltonian.Internal.ExFunctional;
import qchem.ChargeDensity;
 
ExFunctional::ExFunctional()
    : itsChargeDensity(0)
    , isPolarized(true)
{};


void ExFunctional::InsertChargeDensity(const DM_CD* cd)
{
    assert(cd);
    itsChargeDensity=cd;
}

Vector<double> ExFunctional::GetVxcs(const Vector<double>& ros) const
{
    Vector<double> ret(ros.size());
    Vector<double>::iterator i(ret.begin());
    Vector<double>::const_iterator  b(ros.begin());
    for (; i!=ret.end()&&b!=ros.end(); i++,b++) *i=GetVxc(*b);
    return ret;
}


