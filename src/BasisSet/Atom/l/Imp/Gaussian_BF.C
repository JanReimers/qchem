// File: Atom/l/Gaussian_BF.C r^l exp(-a*r^2) type Gaussian basis function.
module;
#include <iostream>
#include <cassert>
#include <cmath>
module qchem.BasisSet.Atom.Internal.l.GaussianBS;
import Common.IntPow;
import oml;

namespace Atoml  
{
namespace Gaussian
{

BasisFunction::BasisFunction()
    : itsExponent     (0)
    , itsL            (0)
    , itsNormalization(0)
{};

BasisFunction::BasisFunction(double theExponent, int theL, double norm)
    : itsExponent     (theExponent)
    , itsL            (theL       )
    , itsNormalization(norm)
{
};

BasisFunction::~BasisFunction() {};

std::ostream& BasisFunction::Write(std::ostream& os) const
{
    return os << itsExponent << " ";
}

double BasisFunction::operator()(const RVec3& r) const
{
    return itsNormalization*uintpow(norm(r),itsL)*exp(-itsExponent*r*r);
}

RVec3 BasisFunction::Gradient(const RVec3& r) const
{
    RVec3 ret(0,0,0);
    double gr=operator()(r);
    double mr=norm(r);
    if (mr>0)
    {
        RVec3 Rhat=r/mr;
        ret = Rhat*gr*(itsL/mr-2.0*mr*itsExponent);
    }
    return ret;
}

}} //namespace