// File: Atom/kappa/Slater_BF.C  Slater basis function with Restricted Kinetic Balance (RKB).

#include "kappa/Slater_BF.H"
#include "Symmetry/Okmj.H"
#include "Common/IntPower.H"
#include "oml/vector3d.h"
#include <iostream>
#include <cassert>

namespace Atom_kappa
{
namespace Slater
{

Large_BasisFunction::Large_BasisFunction(double ex, int _kappa, int _mj, double norm)
    : itsExponent     (ex)
    , kappa(_kappa), mj(_mj), l(Omega_kmj_Sym::l(kappa))
    , itsNormalization(norm)
{
};

std::ostream& Large_BasisFunction::Write(std::ostream& os) const
{
    if (StreamableObject::Pretty())
    {
        os << itsExponent << " ";
    }
    return os;
}

//
//  Y_lm is complex.  To keep it real, just output the radial part for now.
//
double Large_BasisFunction::operator()(const Vec3& r) const
{
    double mr=norm(r);
    return itsNormalization*uintpow(mr,l)*exp(-itsExponent*mr);
}

Large_BasisFunction::Vec3 Large_BasisFunction::Gradient(const Vec3& r) const
{
    Vec3 ret(0,0,0);
    if (r==ret) return ret; //Cusp at the origin so grad is undefined.
    double mr=norm(r);
    assert(mr>0);
    double gr=operator()(r);
    Vec3 r_hat=r/mr;
    return r_hat*gr*(l/mr-itsExponent);
}

::BasisFunction* Large_BasisFunction::Clone() const
{
    return new  Large_BasisFunction(*this);
}

Small_BasisFunction::Small_BasisFunction(const Large_BasisFunction* _Pr,double norm)
: Pr(_Pr)
, itsNormalization(norm)
{
    assert(Pr);
    assert(itsNormalization>0.0);
}

std::ostream& Small_BasisFunction::Write(std::ostream& os) const
{
    return Pr->Write(os);
}

//
//  Y_lm is complex.  To keep it real, just output the radial part for now.
//
double Small_BasisFunction::operator()(const Vec3& r) const
{
    double e=Pr->itsExponent; 
    double f = Pr->kappa >0 ? (2*Pr->kappa+1)/norm(r)-e : -e;
    double n=itsNormalization/Pr->itsNormalization; //Pr(r) is already normalized.
    return n*f*(*Pr)(r); 
}

Large_BasisFunction::Vec3 Small_BasisFunction::Gradient(const Vec3& r) const
{
   assert(false);
    return Pr->Gradient(r);
}

::BasisFunction* Small_BasisFunction::Clone() const
{
    return new  Small_BasisFunction(*this);
}




}} //namespace
