// File: Atom/kappa/Gaussian_IEC.C  Integral Engine Client (IEC) for RKB Gaussians.

#include "Imp/BasisSet/Atom/kappa/Gaussian_IEC.H"
#include "Imp/Symmetry/OkmjQN.H"
#include "Imp/BasisSet/Atom/radial/Gaussian/Integrals.H"

namespace Atom_kappa
{
namespace Gaussian
{

IrrepIEClient::IrrepIEClient(size_t N,int _kappa) 
    : AtomIrrepIEClient(N)
    , kappa(_kappa)
    , j(Omega_kmjQN::j(kappa))
    {};



double IrrepIEClient::Norm(double e, size_t l) const
{
    return ::Gaussian::Norm(e,l);
}

double Small_IrrepIEClient::Norm(double e, size_t l) const
{
    //return GaussianNorm(e,l)/1.0; 
    return 1.0/sqrt(::Gaussian::IE_Primatives::Grad2(e,e,l,l));  //SlaterIntegral already has 4*Pi
}


void Dirac_IrrepIEClient::Init(const IrrepIEClient* liec,const IrrepIEClient* siec)
{
    itsLargeIEC=liec;
    itsSmallIEC=siec;
    assert(itsLargeIEC);
    assert(itsSmallIEC);
}

size_t Dirac_IrrepIEClient::size() const
{
    assert(itsLargeIEC);
    assert(itsSmallIEC);
    return itsLargeIEC->size()+itsSmallIEC->size();
}

}} //namespace
