// File: Atom/l/Gaussian_BF.C r^l exp(-a*r^2) type Gaussian basis function.

#include "Imp/BasisSet/Atom/l/Gaussian_BF.H"
#include "Common/IntPower.H"
#include "oml/vector3d.h"
#include <iostream>
#include <cassert>

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

std::ostream& BasisFunction::Write(std::ostream& os) const
{
    return os << itsExponent << " ";
}

double BasisFunction::operator()(const Vec3& r) const
{
    return itsNormalization*uintpow(norm(r),itsL)*exp(-itsExponent*r*r);
}

BasisFunction::Vec3 BasisFunction::Gradient(const Vec3& r) const
{
    Vec3 ret(0,0,0);
    double gr=operator()(r);
    double mr=norm(r);
    if (mr>0)
    {
        Vec3 Rhat=r/mr;
        ret = Rhat*gr*(itsL/mr-2.0*mr*itsExponent);
    }
    return ret;
}

BasisFunction* BasisFunction::Clone() const
{
    return new  BasisFunction(*this);
}

}} //namespace