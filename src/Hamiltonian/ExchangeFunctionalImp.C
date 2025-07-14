// File: ExchangeFunctionalImplementation.C  Common implementation details for an exchange potential.

#include "oml/vector.h"
#include <cassert>

#include "ExchangeFunctionalImp.H"
#include <ChargeDensity/ChargeDensity.H>

ExFunctionalImp::ExFunctionalImp()
    : itsChargeDensity(0)
    , isPolarized(true)
{};


void ExFunctionalImp::InsertChargeDensity(const DM_CD* cd)
{
    assert(cd);
    itsChargeDensity=cd;
}

void ExFunctionalImp::Eval(const Mesh& m, Vector<double>& v) const
{
    assert(itsChargeDensity);
    Vector<double> ro=(*itsChargeDensity)(m);
    if (!isPolarized) ro*=0.5;
    v=GetVxcs(ro);
}

Vector<double> ExFunctionalImp::GetVxcs(const Vector<double>& ros) const
{
    Vector<double> ret(ros.size());
    Vector<double>::iterator i(ret.begin());
    Vector<double>::const_iterator  b(ros.begin());
    for (; i!=ret.end()&&b!=ros.end(); i++,b++) *i=GetVxc(*b);
    return ret;
}


