
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
    return 1.0/sqrt(2.0*Kinetic(e,e,l));  //SlaterIntegral already has 4*Pi
}

double Small_IrrepIEClient::Kinetic(double ea, double eb,size_t la)
{
    double ab=ea+eb;
    int l=la; //Safer to do formulas with int.
    int ll=l*(l+1);
    double Term1=((l+1)*(l+1)+ll)*SlaterIntegral(ab,2*l-2); //SlaterIntegral already has 4*Pi
    double Term2=-(l+1)*ab* SlaterIntegral(ab,2*l-1);
    double Term3=ea*eb*SlaterIntegral(ab,2*l);
    return 0.5*(Term1+Term2+Term3);
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
