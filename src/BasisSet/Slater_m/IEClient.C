
#include "Imp/BasisSet/Slater_m/IEClient.H"
#include "Imp/Integrals/SlaterCD.H"

namespace Slater_m
{
    

const Cacheable* IEClient::Create(size_t ia,size_t ic,size_t ib,size_t id) const
{
    return new SlaterCD(unique_esv[ia]+unique_esv[ib],unique_esv[ic]+unique_esv[id],LMax(ia,ib,ic,id));
}

Vector<double> IEClient::loop_4_direct(size_t id, size_t la, size_t lc)  const
{
    const Cacheable* c=Cache4::loop_4(id);
    const SlaterCD* cd = dynamic_cast<const SlaterCD*>(c);
    return cd->Coulomb_Rk(la,lc);
}
Vector<double> IEClient::loop_4_exchange(size_t id, size_t la, size_t lc)  const
{
    const Cacheable* c=Cache4::loop_4(id);
    const SlaterCD* cd = dynamic_cast<const SlaterCD*>(c);
    return cd->ExchangeRk(la,lc);
}
    
} //namespace
