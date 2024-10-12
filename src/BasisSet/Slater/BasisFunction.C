// File: Slater.C  r^(n-1) exp(-ar) type basis function for an atom.

#include "Imp/BasisSet/Slater/BasisFunction.H"
#include "Misc/IntPower.H"
#include "oml/vector3d.h"
#include "oml/imp/binio.h"
#include <cmath>
#include <iostream>
#include <cassert>

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
    UniqueID::Write(os);
    if ( StreamableObject::Binary())
    {
        BinaryWrite(itsExponent     ,os);
        BinaryWrite(itsN            ,os);
        BinaryWrite(itsL            ,os);
        BinaryWrite(itsNormalization,os);
    }
    if (StreamableObject::Ascii ()) os << itsExponent << " " << itsN << " " << itsL << " " << itsNormalization << " ";
    if (StreamableObject::Pretty())
    {
//    if (itsL >0) os << "r"; else os << " ";
//    if (itsL >1) os << "^" << itsL; else os << "  ";
//    os << " exp(-" << itsExponent << "*r)" << std::endl;
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
        BinaryRead(itsN            ,is);
        BinaryRead(itsL            ,is);
        BinaryRead(itsNormalization,is);
    }
    else
    {
        is >> itsExponent >> itsN >> itsL >> itsNormalization ;
        is.get();
    }
    return is;
}

double BasisFunction::operator()(const Vec3& r) const
{
    double mr=norm(r);
    return itsNormalization*uintpow(mr,itsN-1)*exp(-itsExponent*mr);
}

BasisFunction::Vec3 BasisFunction::Gradient(const Vec3& r) const
{
    assert(r!=Vec3(0,0,0)); //cusp at the origin.
    double mr=norm(r);
    assert(mr>0);
    Vec3 r_hat=r/mr;
    double g=(*this)(r);
    return r_hat*g*((itsN-1)/mr-itsExponent);
}

BasisFunction* BasisFunction::Clone() const
{
    return new  BasisFunction(*this);
}

} //namespace
