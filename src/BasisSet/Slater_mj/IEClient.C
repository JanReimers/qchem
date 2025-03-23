
#include "Imp/BasisSet/Slater_mj/IEClient.H"
#include "Imp/Symmetry/OkmjQN.H"
#include "Imp/Integrals/SlaterIntegrals.H"
using std::cout;
using std::endl;

namespace Slater_mj
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
     return SlaterNorm(e,l+1);  //Already has 1/sqrt(4*Pi).
}

double Small_IrrepIEClient::Norm(double e, size_t l) const
{ 
    return 1.0/sqrt(Slater::IE_Primatives::Grad2(e,e,l,l));  //SlaterIntegral already has 4*Pi
}
    
void Dirac_IrrepIEClient::Init(const Slater_mj::IrrepIEClient* liec,const Slater_mj::IrrepIEClient* siec)
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
