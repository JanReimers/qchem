
#include "Imp/BasisSet/Slater/IEClient.H"
#include "Imp/BasisSet/Atom/radial/Slater/Integrals.H"



namespace Slater
{
    
double IrrepIEClient::Norm(double e, size_t l) const
{
     return SlaterNorm(e,l+1);  
}





} //namespace
