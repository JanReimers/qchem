
#include "Imp/BasisSet/Slater_mj/IEClient.H"
#include "Imp/Symmetry/OkmjQN.H"
#include "Imp/Integrals/SlaterIntegrals.H"


namespace Slater_mj
{
IrrepIEClient::IrrepIEClient(size_t N,int _kappa, double _mj) 
    : AtomIrrepIEClient(N)
    , kappa(_kappa)
    , j(Omega_kmjQN::j(kappa))
    , mj(_mj) 
    {};
    
void IrrepIEClient::Init(const Vector  <double>& exponents)
{
    AtomIrrepIEClient::Init(exponents,Omega_kmjQN::l(kappa),Omega_kmjQN::ml(kappa,mj));
}

double IrrepIEClient::Norm(double e, size_t l) const
{
     return SlaterNorm(e,l+1);  
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

} //namespace
