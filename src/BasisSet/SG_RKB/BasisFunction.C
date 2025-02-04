// File: SphericalGaussian_RKB/BF.C  r^l exp(-ar^2) type basis function for an atom.

#include "Imp/BasisSet/SG_RKB/BasisFunction.H"
#include "Imp/Symmetry/OkmjQN.H"
#include "Imp/Misc/IntPower.H"
#include "oml/vector3d.h"
#include "oml/imp/binio.h"
#include <cmath>
#include <iostream>
#include <cassert>

namespace SphericalGaussian_RKB
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
    , l               (Omega_kmjQN::l(kappa))       
    , itsNormalization(norm)
{
};

bool Large_BasisFunction::operator==(const ::BasisFunction& bf) const
{
    const Large_BasisFunction& sbf = dynamic_cast<const Large_BasisFunction&>(bf);
    assert(&sbf);
    return itsExponent==(sbf.itsExponent) && (kappa==sbf.kappa) && (l==sbf.l);
}

std::ostream& Large_BasisFunction::Write(std::ostream& os) const
{
    if (StreamableObject::Pretty())
    {
        os << itsExponent << " ";
    }
    return os;
}

std::istream& Large_BasisFunction::Read(std::istream& is)
{
    return is;
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

bool Small_BasisFunction::operator==(const ::BasisFunction& bf) const
{
    return *Pr==bf;
}

std::ostream& Small_BasisFunction::Write(std::ostream& os) const
{
    return Pr->Write(os);
}

std::istream& Small_BasisFunction::Read(std::istream& is)
{
    return is;
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


} //namespace
