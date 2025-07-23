// File: Atom/kappa/Gaussian_BF.C  Gaussians with Restricted Kinetic Balance (RKB).
module;
#include <cmath>
#include <iostream>
#include <cassert>
module qchem.BasisSet.Atom.kappa.GaussianBS;
import Common.IntPow;
import qchem.Symmetry.Okmj;
import oml;

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
    os << itsExponent << " ";
    return os;
}

double Large_BasisFunction::operator()(const RVec3& r) const
{
    return itsNormalization*uintpow(norm(r),l)*exp(-itsExponent*r*r);
}

RVec3 Large_BasisFunction::Gradient(const RVec3& r) const
{
    RVec3 ret(0,0,0);
    double gr=operator()(r);
    double mr=norm(r);
    if (mr>0)
    {
        RVec3 Rhat=r/mr;
        ret = Rhat*gr*(l/mr-2.0*mr*itsExponent);
    }
    return ret;
}

Real_BF* Large_BasisFunction::Clone() const
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
double Small_BasisFunction::operator()(const RVec3& r) const
{
    double e=Pr->itsExponent; 
    double mr=norm(r);
    double f = Pr->kappa >0 ? (2*Pr->kappa+1)/mr-2*e*mr : -2*e*mr;
    double n=itsNormalization/Pr->itsNormalization; //Pr(r) is already normalized.
    return n*f*(*Pr)(r); 
}

RVec3 Small_BasisFunction::Gradient(const RVec3& r) const
{
   assert(false);
    return Pr->Gradient(r);
}

::Real_BF* Small_BasisFunction::Clone() const
{
    return new  Small_BasisFunction(*this);
}


}} //namespace
