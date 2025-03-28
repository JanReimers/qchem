// File: Atom/l/Slater_IEC.C Integral Engine Client for Slater atom basis functions.

#include "Imp/BasisSet/Atom/l/Slater_IEC.H"
#include "Imp/BasisSet/Atom/radial/Slater/Integrals.H"


namespace Atoml
{
namespace Slater
{
    
double IrrepIEClient::Norm(double e, size_t l) const
{
     return ::Slater::Norm(e,l+1);  
}





}} //namespace
