
#include "Imp/BasisSet/Slater_m/IEClient.H"

namespace Slater_m
{
    

const Cacheable* IEClient::Create(size_t ia,size_t ic,size_t ib,size_t id) const
{
    return new SlaterCD(unique_esv[ia]+unique_esv[ib],unique_esv[ic]+unique_esv[id],LMax());
}

} //namespace
