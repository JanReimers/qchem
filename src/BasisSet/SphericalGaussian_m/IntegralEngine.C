// File: SphericalGaussianIE1.C  Here is where all the integral get calculated.


#include "Imp/BasisSet/SphericalGaussian_m/IntegralEngine.H"
#include "Imp/Integrals/AngularIntegrals.H"

namespace SphericalGaussian_m
{
    
IntegralEngine::RVec IntegralEngine::Coulomb_AngularIntegrals(size_t la, size_t lc, int ma, int mc) const
{
    return AngularIntegrals::Coulomb(la,lc,ma,mc);
}

IntegralEngine::RVec IntegralEngine::ExchangeAngularIntegrals(size_t la, size_t lb, int ma, int mb) const
{
    return AngularIntegrals::Exchange(la,lb,ma,mb);
}


} //namespace
