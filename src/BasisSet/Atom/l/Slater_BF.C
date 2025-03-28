// File: Atom/l/Slater_BF.H  r^l exp(-a*r) type Slater basis function.

#include "Imp/BasisSet/Atom/l/Slater_BF.H"
#include "Imp/Misc/IntPower.H"
#include "oml/vector3d.h"
#include <iostream>
#include <cassert>

namespace Atoml
{
namespace Slater
{

BasisFunction::BasisFunction()
    : itsExponent     (0)
    , itsN            (0)
    , itsL            (0)
    , itsNormalization(0)
{};

BasisFunction::BasisFunction(double ex, int n, int l, double norm)
    : itsExponent     (ex)
    , itsN            (n)
    , itsL            (l)
    , itsNormalization(norm)
{
};

bool BasisFunction::operator==(const ::BasisFunction& bf) const
{
    const BasisFunction& sbf = dynamic_cast<const BasisFunction&>(bf);
    assert(&sbf);
    return itsExponent==(sbf.itsExponent) && (itsN==sbf.itsN) && (itsL==sbf.itsL);
}

std::ostream& BasisFunction::Write(std::ostream& os) const
{
    return os << itsExponent << " ";
}

double BasisFunction::operator()(const Vec3& r) const
{
    double mr=norm(r);
    return itsNormalization*uintpow(mr,itsN-1)*exp(-itsExponent*mr);
}

BasisFunction::Vec3 BasisFunction::Gradient(const Vec3& r) const
{
    Vec3 ret(0,0,0);
    if (r==ret) return ret; //Cusp at the origin so grad is undefined.
    double mr=norm(r);
    assert(mr>0);
    assert(itsL==itsN-1);
    double gr=operator()(r);
    Vec3 r_hat=r/mr;
    return r_hat*gr*(itsL/mr-itsExponent);
}

BasisFunction* BasisFunction::Clone() const
{
    return new  BasisFunction(*this);
}

}} //namespace
