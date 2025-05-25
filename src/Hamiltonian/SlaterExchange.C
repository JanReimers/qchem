// File: SlaterExchange.C  Slater exchange potential.



#include "Imp/Hamiltonian/SlaterExchange.H"
#include <ChargeDensity.H>
#include "oml/imp/binio.h"
#include <iostream>
#include <cmath>
#include <cassert>

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

double SlaterExchange::operator()(const Vec3& r) const
{
    double ro = (*itsChargeDensity)(r);
    return GetVxc(ro);
}

double SlaterExchange::GetVxc(double ro) const
{
    if (itsSpin==Spin::None) ro*=0.5;
    double ret=0;
    if (ro > 0.0)
    {
        ret=-3.0 * itsAlpha * pow(3.0*ro/(M_PI*4.0) , 1.0/3.0);
    }
    return ret;
}

SlaterExchange::Vec3 SlaterExchange::Gradient(const Vec3& r) const
{
    double ro = (*itsChargeDensity)(r);
    if (itsSpin==Spin::None) ro*=0.5;
    Vec3 ret(0,0,0);
    if (ro > 0.0)
    {
        ret=-3.0 * itsAlpha * pow(3.0*ro/(M_PI*4.0) , -2.0/3.0) / (M_PI*4.0) * itsChargeDensity->Gradient(r);
    }
    if (!isPolarized) ret*=0.5;
    return ret;
}

std::ostream& SlaterExchange::Write(std::ostream& os) const
{
    if (StreamableObject::Binary())
    {
        BinaryWrite(itsAlpha,os);
    }
    else
    {
        os << itsAlpha << " ";
    }
    return os;
}

std::istream& SlaterExchange::Read (std::istream& is)
{
    if (StreamableObject::Binary())
    {
        BinaryRead(itsAlpha,is);
    }
    else
    {
        is >> itsAlpha;
        assert(is.get() == ' ');
    }
    return is;
}


