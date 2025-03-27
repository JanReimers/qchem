// File: Atom/ml/Gaussian_BF.C  r^l exp(-ar^2)*Y_lm type basis function.

#include "Imp/BasisSet/Atom/ml/Gaussian_BF.H"
#include "Imp/Misc/IntPower.H"
#include "oml/vector3d.h"
// #include <cmath>
#include <iostream>
#include <cassert>

namespace Atom_ml
{
namespace Gaussian
{
    
//#######################################################################
//
//  Spherical gaussian implementation, adds orbital AM Q#.
//

BasisFunction::BasisFunction()
    : itsExponent     (0)
    , n(0), l(0), m(0)
    , itsNormalization(0)
{};

BasisFunction::BasisFunction(double ex, int _n, int _l, int _m, double norm)
    : itsExponent     (ex)
    , n(_n), l(_l), m(_m)
    , itsNormalization(norm)
{
};

bool BasisFunction::operator==(const ::BasisFunction& bf) const
{
    const BasisFunction& sgbf = dynamic_cast<const BasisFunction&>(bf);
    assert(&sgbf);
    return itsExponent==(sgbf.itsExponent) && (n==sgbf.n) && (l==sgbf.l) && (m==sgbf.m);
}

std::ostream& BasisFunction::Write(std::ostream& os) const
{
    return os << itsExponent << " ";
}

std::istream& BasisFunction::Read(std::istream& is)
{
    return is;
}

double BasisFunction::operator()(const Vec3& r) const
{
    return itsNormalization*uintpow(norm(r),n-1)*exp(-itsExponent*r*r);
}

BasisFunction::Vec3 BasisFunction::Gradient(const Vec3& r) const
{
    Vec3 ret(0,0,0);
    double gr=operator()(r);
    double mr=norm(r);
    if (mr>0)
    {
        Vec3 Rhat=r/mr;
        ret = Rhat*gr*((n-1)/mr-2.0*mr*itsExponent);
    }
    return ret;
}

BasisFunction* BasisFunction::Clone() const
{
    return new  BasisFunction(*this);
}

} //namespace
} //namespace
