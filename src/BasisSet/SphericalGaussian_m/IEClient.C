
#include "Imp/BasisSet/SphericalGaussian_m/IEClient.H"
#include "Imp/Integrals/SphericalGaussianCD.H"

namespace SphericalGaussian_m
{


const Cacheable* IEClient::Create(size_t ia,size_t ic,size_t ib,size_t id) const
{
    return new SphericalGaussianCD(unique_esv[ia]+unique_esv[ib],unique_esv[ic]+unique_esv[id],LMax(ia,ib,ic,id));
}

Vector<double> IEClient::loop_4_direct(size_t id, size_t la, size_t lc)  const
{
    const Cacheable* c=Cache4::loop_4(es_indices[id-1]);
    const SphericalGaussianCD* cd = dynamic_cast<const SphericalGaussianCD*>(c);
    return cd->Coulomb_Rk(la,lc);
}
Vector<double> IEClient::loop_4_exchange(size_t id, size_t la, size_t lc)  const
{
    const Cacheable* c=Cache4::loop_4(es_indices[id-1]);
    const SphericalGaussianCD* cd = dynamic_cast<const SphericalGaussianCD*>(c);
    return cd->ExchangeRk(la,lc);
}


} //namespace
