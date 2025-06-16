// File: Atom/kappa/Gaussian_BF.C  Gaussians with Restricted Kinetic Balance (RKB).

#include "Atom/kappa/Gaussian_BF.H"
#include "Symmetry/Okmj.H"
#include "Common/IntPower.H"
#include "oml/vector3d.h"
#include "oml/imp/binio.h"
#include "oml/imp/stream.h"
#include <cmath>
#include <iostream>
#include <cassert>

namespace Atom_kappa
{
namespace Gaussian
{
//#######################################################################
//
//  Spherical gaussian implementation, adds orbital AM Q#.
//

Large_BasisFunction::Large_BasisFunction()
    : itsExponent     (0)
    , kappa           (0)
    , l               (0)
    , itsNormalization(0)
{};

Large_BasisFunction::Large_BasisFunction(double theExponent, int _kappa, double norm)
    : itsExponent     (theExponent)
    , kappa           (_kappa)
    , l               (Omega_kmj_Sym::l(kappa))       
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

double Large_BasisFunction::operator()(const Vec3& r) const
{
    return itsNormalization*uintpow(norm(r),l)*exp(-itsExponent*r*r);
}

Large_BasisFunction::Vec3 Large_BasisFunction::Gradient(const Vec3& r) const
{
    Vec3 ret(0,0,0);
    double gr=operator()(r);
    double mr=norm(r);
    if (mr>0)
    {
        Vec3 Rhat=r/mr;
        ret = Rhat*gr*(l/mr-2.0*mr*itsExponent);
    }
    return ret;
}

BasisFunction* Large_BasisFunction::Clone() const
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
    double mr=norm(r);
    double f = Pr->kappa >0 ? (2*Pr->kappa+1)/mr-2*e*mr : -2*e*mr;
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
