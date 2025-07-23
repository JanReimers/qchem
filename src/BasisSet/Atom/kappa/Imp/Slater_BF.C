// File: Atom/kappa/Slater_BF.C  Slater basis function with Restricted Kinetic Balance (RKB).
module;
#include <iostream>
#include <cassert>
#include <cmath>
module qchem.BasisSet.Atom.kappa.SlaterBS;
import Common.IntPow;
import qchem.Symmetry.Okmj;
import oml;

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
    os << itsExponent << " ";
    return os;
}

//
//  Y_lm is complex.  To keep it real, just output the radial part for now.
//
double Large_BasisFunction::operator()(const RVec3& r) const
{
    double mr=norm(r);
    return itsNormalization*uintpow(mr,l)*exp(-itsExponent*mr);
}

RVec3 Large_BasisFunction::Gradient(const RVec3& r) const
{
    RVec3 ret(0,0,0);
    if (r==ret) return ret; //Cusp at the origin so grad is undefined.
    double mr=norm(r);
    assert(mr>0);
    double gr=operator()(r);
    RVec3 r_hat=r/mr;
    return r_hat*gr*(l/mr-itsExponent);
}

::Real_BF* Large_BasisFunction::Clone() const
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
    double f = Pr->kappa >0 ? (2*Pr->kappa+1)/norm(r)-e : -e;
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
