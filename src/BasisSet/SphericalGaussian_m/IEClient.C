
#include "Imp/BasisSet/SphericalGaussian_m/IEClient.H"

namespace SphericalGaussian_m
{


const Cacheable* IEClient::Create(size_t ia,size_t ic,size_t ib,size_t id) const
{
    return new SphericalGaussianCD(unique_esv[ia]+unique_esv[ib],unique_esv[ic]+unique_esv[id],LMax());
}




} //namespace
