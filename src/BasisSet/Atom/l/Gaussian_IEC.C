// File: Atom/l/Gaussian_IEC.H

#include "Imp/BasisSet/Atom/l/Gaussian_IEC.H"
#include "Imp/BasisSet/Atom/radial/Gaussian/Integrals.H"

namespace Atoml::Gaussian
{

double IrrepIEClient::Norm(double e, size_t l) const
{
    return ::Gaussian::Norm(e,l);
}


} //namespace
