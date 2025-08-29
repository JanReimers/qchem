// File: BasisSet/Atom/radial/Imp/Gaussian_BS.C
module;
#include <cassert>
module BasisSet.Atom.Gaussian_BS;
import qchem.BasisSet.Atom.Internal.radial.Gaussian.Rk;

void Gaussian_BS::Register(IBS_Evaluator * eval)
{
    assert(eval);
    eval->Register(&grouper);
}

Rk* Gaussian_BS::Create(size_t ia,size_t ic,size_t ib,size_t id) const
{
    return new Gaussian::RkEngine(
        grouper.unique_esv[ia]+grouper.unique_esv[ib],
        grouper.unique_esv[ic]+grouper.unique_esv[id],
        grouper.LMax(ia,ib,ic,id));
}

RVec Gaussian_BS::loop_4_direct(size_t id, size_t la, size_t lc)  const
{
    const Cacheable* c=Cache4::loop_4(id);
    const Gaussian::RkEngine* cd = dynamic_cast<const Gaussian::RkEngine*>(c);
    return cd->Coulomb_Rk(la,lc);
}
RVec Gaussian_BS::loop_4_exchange(size_t id, size_t la, size_t lc)  const
{
    const Cacheable* c=Cache4::loop_4(id);
    const Gaussian::RkEngine* cd = dynamic_cast<const Gaussian::RkEngine*>(c);
    return cd->ExchangeRk(la,lc);
}
    