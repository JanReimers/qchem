
#include "Imp/BasisSet/SG_RKB/IEClient.H"
#include "Imp/Symmetry/OkmjQN.H"
#include "Imp/Integrals/GaussianIntegrals.H"
#include "Imp/Integrals/SphericalGaussianCD.H"

namespace SphericalGaussian_RKB
{

IrrepIEClient::IrrepIEClient(size_t N,int _kappa) 
    : AtomIrrepIEClient(N)
    , kappa(_kappa)
    , j(Omega_kmjQN::j(kappa))
    {};

void IrrepIEClient::Init(const Vector  <double>& exponents)
{
    AtomIrrepIEClient::Init(exponents,Omega_kmjQN::l(kappa));
}

double IrrepIEClient::Norm(double e, size_t l) const
{
    return GaussianNorm(e,l);
}

double Small_IrrepIEClient::Norm(double e, size_t l) const
{
    //return GaussianNorm(e,l)/1.0; 
    return 1.0/sqrt(2.0*SphericalGaussian::IE_Primatives::Kinetic(e,e,l,l));  //SlaterIntegral already has 4*Pi
}


void Dirac_IrrepIEClient::Init(const SphericalGaussian_RKB::IrrepIEClient* liec,const SphericalGaussian_RKB::IrrepIEClient* siec)
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

} //namespace
