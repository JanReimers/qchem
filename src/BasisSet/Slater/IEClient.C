
#include "Imp/BasisSet/Slater/IEClient.H"
#include "Imp/Integrals/SlaterIntegrals.H"
#include "Imp/Integrals/SlaterCD.H"



namespace Slater
{
    
double IrrepIEClient::Norm(double e, size_t l) const
{
     return SlaterNorm(e,l+1);  
}





} //namespace
