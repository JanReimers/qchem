// File: SphericalGaussian.C  r^l exp(-ar^2) type basis function for an atom.



#include "Imp/BasisSet/Atom/radial/Gaussian/BasisFunction.H"
#include "Imp/Misc/IntPower.H"
#include "oml/vector3d.h"
#include <iostream>
#include <cassert>

namespace Gaussian
{
//#######################################################################
//
//  Spherical gaussian function implementation.
//

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

bool BasisFunction::operator==(const ::BasisFunction& bf) const
{
    const BasisFunction& sgbf = dynamic_cast<const BasisFunction&>(bf);
    assert(&sgbf);
    return itsExponent==(sgbf.itsExponent) && (itsL==sgbf.itsL);
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

BasisFunction_ml::BasisFunction_ml(double e,int n, int l, int _ml, double norm)
: BasisFunction(e,l,norm)
, ml(_ml)
{};

bool BasisFunction_ml::operator==(const ::BasisFunction& bf) const
{
    const BasisFunction_ml& sgbf = dynamic_cast<const BasisFunction_ml&>(bf);
    assert(&sgbf);
    return  BasisFunction::operator==(bf) && (ml==sgbf.ml);
}

BasisFunction* BasisFunction_ml::Clone() const
{
    return new  BasisFunction_ml(*this);
}

} //namespace
