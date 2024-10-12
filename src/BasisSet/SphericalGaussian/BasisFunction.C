// File: SphericalGaussian.C  r^l exp(-ar^2) type basis function for an atom.



#include "Imp/BasisSet/SphericalGaussian/BasisFunction.H"
#include "Imp/Integrals/GaussianIntegrals.H"
#include "Misc/IntPower.H"
#include "oml/vector3d.h"
#include "oml/imp/binio.h"
#include <cmath>
#include <iostream>
#include <cassert>

namespace SphericalGaussian
{
//#######################################################################
//
//  Spherical gaussian implementation, adds orbital AM Q#.
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
    UniqueID::Write(os);
    if ( StreamableObject::Binary())
    {
        BinaryWrite(itsExponent     ,os);
        BinaryWrite(itsL            ,os);
        BinaryWrite(itsNormalization,os);
    }
    if (StreamableObject::Ascii ()) os << itsExponent << " " << itsL << " " << itsNormalization << " ";
    if (StreamableObject::Pretty())
    {
//    if (itsL >0) os << "r"; else os << " ";
//    if (itsL >1) os << "^" << itsL; else os << "  ";
//    os << " exp(-" << itsExponent << " r^2)" << std::endl;
        os << itsExponent << " ";
    }
    return os;
}

std::istream& BasisFunction::Read(std::istream& is)
{
    UniqueID::Read(is);
    if (StreamableObject::Binary())
    {
        BinaryRead(itsExponent     ,is);
        BinaryRead(itsL            ,is);
        BinaryRead(itsNormalization,is);
    }
    else
    {
        is >> itsExponent >> itsL >> itsNormalization ;
        is.get();
    }
    return is;
}

double BasisFunction::operator()(const Vec3& r) const
{
    return itsNormalization*uintpow(norm(r),itsL)*exp(-itsExponent*r*r);
}

BasisFunction::Vec3 BasisFunction::Gradient(const Vec3& r) const
{
    Vec3 ret(0,0,0);
    double mr=norm(r);
    if (mr>0)
    {
        double diff = itsL==0 ? 0 : uintpow(norm(r),itsL-1);
        ret = -2.0*r*itsExponent*(*this)(r) + itsNormalization*normalize(r)*diff*exp(-itsExponent*r*r);
    }
    return ret;
}

BasisFunction* BasisFunction::Clone() const
{
    return new  BasisFunction(*this);
}

} //namespace
