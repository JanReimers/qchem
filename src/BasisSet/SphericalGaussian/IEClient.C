
#include "Imp/BasisSet/SphericalGaussian/IEClient.H"
#include "Imp/Integrals/GaussianIntegrals.H"
#include "Imp/Integrals/SphericalGaussianCD.H"

namespace SphericalGaussian
{

double IrrepIEClient::Norm(double e, size_t l) const
{
    return GaussianNorm(e,l);
}


} //namespace
