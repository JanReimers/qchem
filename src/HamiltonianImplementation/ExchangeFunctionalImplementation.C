// File: ExchangeFunctionalImplementation.C  Common implementation details for an exchange potential.



#include "HamiltonianImplementation/ExchangeFunctionalImplementation.H"
#include "ChargeDensity.H"
#include <cassert>

ExchangeFunctionalImplementation::ExchangeFunctionalImplementation()
    : ScalarFunctionBuffer<double>(false,false)
    , itsChargeDensity(0)
    , isPolarized(true)
{};


void ExchangeFunctionalImplementation::InsertChargeDensity(const ChargeDensity* theChargeDensity)
{
    assert(theChargeDensity);
    itsChargeDensity=theChargeDensity;
    MakeBufferDirty();
}

void ExchangeFunctionalImplementation::Eval(const Mesh& m, Vector<double>& v) const
{
    assert(itsChargeDensity);
    Vector<double> ro=(*itsChargeDensity)(m);
    if (!isPolarized) ro*=0.5;
    v=GetVxcs(ro);
}

Vector<double> ExchangeFunctionalImplementation::GetVxcs(const Vector<double>& ros) const
{
    Vector<double> ret(ros.size());
    Vector<double>::iterator i(ret.begin());
    Vector<double>::const_iterator  b(ros.begin());
    for (; i!=ret.end()&&b!=ros.end(); i++,b++) *i=GetVxc(*b);
    return ret;
}


