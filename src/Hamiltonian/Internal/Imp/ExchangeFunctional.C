// File: Imp/ExchangeFunctional.C   Exchange potential for DFT.
module;
#include <cassert>
#include "blaze/Math.h"
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

rvec_t ExFunctional::GetVxcs(const Vector<double>& ros) const
{
    rvec_t ret(ros.size());
    auto i(ret.begin());
    auto b(ros.begin());
    for (; i!=ret.end()&&b!=ros.end(); i++,b++) *i=GetVxc(*b);
    return ret;
}


