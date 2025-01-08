
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

double Small_IrrepIEClient::Kinetic(double ea, double eb,size_t l)
{
    double ab=ea+eb;
    int na=l+1,nb=l+1;
    int n=na+nb;
    double Term1=0.5*(na*nb+l*(l+1))*SlaterIntegral(ab,n-2); //SlaterIntegral already has 4*Pi
    double Term2=-0.5*(na*eb+nb*ea)* SlaterIntegral(ab,n-1);
    double Term3=0.5*ea*eb*SlaterIntegral(ab,n);
    //cout << "Small_IrrepIEClient::Kinetic Terms 1,2,3=" << Term1 << " " << Term2 << " " << Term3 << endl;
    return (Term1+Term2+Term3);
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
