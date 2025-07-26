// File: Atom/l/Slater_BF.H  r^l exp(-a*r) type Slater basis function.
module;
#include <iostream>
#include <cassert>
#include <cmath>
module qchem.BasisSet.Atom.Internal.l.SlaterBS;
import Common.IntPow;
import oml;

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

std::ostream& BasisFunction::Write(std::ostream& os) const
{
    return os << itsExponent << " ";
}

double BasisFunction::operator()(const RVec3& r) const
{
    double mr=norm(r);
    return itsNormalization*uintpow(mr,itsN-1)*exp(-itsExponent*mr);
}

RVec3 BasisFunction::Gradient(const RVec3& r) const
{
    RVec3 ret(0,0,0);
    if (r==ret) return ret; //Cusp at the origin so grad is undefined.
    double mr=norm(r);
    assert(mr>0);
    assert(itsL==itsN-1);
    double gr=operator()(r);
    RVec3 r_hat=r/mr;
    return r_hat*gr*(itsL/mr-itsExponent);
}

BasisFunction* BasisFunction::Clone() const
{
    return new  BasisFunction(*this);
}

}} //namespace
