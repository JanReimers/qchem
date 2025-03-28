// File: Atom/kappa/Slater_IEC.C  Inegral Engine Client (IEC) for Slater basis set with Restricted Kinetic Balance (RKB).

#include "Imp/BasisSet/Atom/kappa/Slater_IEC.H"
#include "Imp/Symmetry/OkmjQN.H"
#include "Imp/BasisSet/Atom/radial/Slater/Integrals.H"
using std::cout;
using std::endl;

namespace Atom_kappa
{
namespace Slater
{

IrrepIEClient::IrrepIEClient(size_t N,int _kappa) 
    : AtomIrrepIEClient(N)
    , kappa(_kappa)
    , j(Omega_kmjQN::j(kappa))
    {};
    


double IrrepIEClient::Norm(double e, size_t l) const
{
     return ::Slater::Norm(e,l+1);  //Already has 1/sqrt(4*Pi).
}

double Small_IrrepIEClient::Norm(double e, size_t l) const
{ 
    return 1.0/sqrt(::Slater::IE_Primatives::Grad2(e,e,l,l));  //SlaterIntegral already has 4*Pi
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
