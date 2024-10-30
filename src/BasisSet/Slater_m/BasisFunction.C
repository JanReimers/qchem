// File: Slater_m/BasisFunction.C  r^l exp(-ar^2) type basis function for an Spherical, with m QN.

#include "Imp/BasisSet/Slater_m/BasisFunction.H"
#include "Imp/Misc/IntPower.H"
#include "oml/vector3d.h"
#include "oml/imp/binio.h"
#include <cmath>
#include <iostream>
#include <cassert>

namespace Slater_m
{

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
    const BasisFunction& sbf = dynamic_cast<const BasisFunction&>(bf);
    assert(&sbf);
    return itsExponent==(sbf.itsExponent) && (n==sbf.n) && (l==sbf.l) && (m==sbf.m);
}

std::ostream& BasisFunction::Write(std::ostream& os) const
{
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
    return is;
}

//
//  Y_lm is complex.  To keep it real, just output the radial part for now.
//
double BasisFunction::operator()(const Vec3& r) const
{
    double mr=norm(r);
    return itsNormalization*uintpow(mr,n-1)*exp(-itsExponent*mr);
}

BasisFunction::Vec3 BasisFunction::Gradient(const Vec3& r) const
{
    Vec3 ret(0,0,0);
    if (r==ret) return ret; //Cusp at the origin so grad is undefined.
    double mr=norm(r);
    assert(mr>0);
    assert(l==n-1);
    double gr=operator()(r);
    Vec3 r_hat=r/mr;
    return r_hat*gr*(l/mr-itsExponent);
}

BasisFunction* BasisFunction::Clone() const
{
    return new  BasisFunction(*this);
}

} //namespace
