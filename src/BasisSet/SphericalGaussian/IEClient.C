
#include "Imp/BasisSet/SphericalGaussian/IEClient.H"
#include "Imp/BasisSet/Atom/radial/Gaussian/Integrals.H"

namespace SphericalGaussian
{

double IrrepIEClient::Norm(double e, size_t l) const
{
    return GaussianNorm(e,l);
}


} //namespace
