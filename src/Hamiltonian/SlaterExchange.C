// File: SlaterExchange.C  Slater exchange potential.



#include "SlaterExchange.H"
#include <ChargeDensity/ChargeDensity.H>
#include <iostream>
#include <cmath>
#include <cassert>

import Common.Constants;


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
        ret=-3.0 * itsAlpha * pow(3.0*ro/FourPi , 1.0/3.0);
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
        ret=-3.0 * itsAlpha * pow(3.0*ro/FourPi , -2.0/3.0) / FourPi * itsChargeDensity->Gradient(r);
    }
    if (!isPolarized) ret*=0.5;
    return ret;
}

std::ostream& SlaterExchange::Write(std::ostream& os) const
{
    os << itsAlpha << " ";
    return os;
}


