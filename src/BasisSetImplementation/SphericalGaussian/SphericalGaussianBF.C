// File: SphericalGaussian.C  r^l exp(-ar^2) type basis function for an atom.



#include "BasisSetImplementation/SphericalGaussian/SphericalGaussianBF.H"
#include "BasisSetImplementation/SphericalGaussian/GaussianIntegrals.H"
#include "Misc/IntPower.H"
#include "oml/vector3d.h"
#include "oml/imp/binio.h"
#include <cmath>
#include <iostream>
#include <cassert>

//#######################################################################
//
//  Spherical gaussian implementation, adds orbital AM Q#.
//

SphericalGaussianBF::SphericalGaussianBF()
    : itsExponent     (0)
    , itsL            (0)
    , itsNormalization(0)
    , itsCharge       (0)
{};

SphericalGaussianBF::SphericalGaussianBF(double theExponent, int theL)
    : itsExponent     (theExponent)
    , itsL            (theL       )
    , itsNormalization(0)
    , itsCharge       (0)
{
};

void SphericalGaussianBF::Init(double norm, double charge)
{   
    assert(!std::isnan(norm));
    assert(norm>0);
    assert(charge>0);
//  Ascii IO messes this up.
//    std::cout << norm << " " << GaussianNorm(itsExponent,itsL) << std::endl;
//    std::cout << log10(fabs(norm-GaussianNorm(itsExponent,itsL))) << std::endl;
//    assert(fabs(norm-GaussianNorm(itsExponent,itsL))<1e-13);
//    assert(fabs(charge-norm*GaussianIntegral(itsExponent,itsL))<1e-14);

    itsNormalization=norm;
    itsCharge=charge;
}

bool SphericalGaussianBF::operator==(const BasisFunction& bf) const
{
    const SphericalGaussianBF& sgbf = dynamic_cast<const SphericalGaussianBF&>(bf);
    assert(&sgbf);
    return itsExponent==(sgbf.itsExponent) && (itsL==sgbf.itsL);
}

double SphericalGaussianBF::GetNormalization() const
{
    return itsNormalization;
}

double SphericalGaussianBF::GetCharge() const
{
    return itsCharge;
}

std::ostream& SphericalGaussianBF::Write(std::ostream& os) const
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

std::istream& SphericalGaussianBF::Read(std::istream& is)
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

double SphericalGaussianBF::operator()(const Vec3& r) const
{
    double a=TwoLPlusOne(itsL);
    return itsNormalization*sqrt(a)*uintpow(norm(r),itsL)*exp(-itsExponent*r*r);
}

SphericalGaussianBF::Vec3 SphericalGaussianBF::Gradient(const Vec3& r) const
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

BasisFunction* SphericalGaussianBF::Clone() const
{
    return new  SphericalGaussianBF(*this);
}

