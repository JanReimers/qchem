// File: BasisSet/Atom/radial/Imp/Slater_BS.C
module;
#include <cassert>
module BasisSet.Atom.Slater.NR.BS_Evaluator;
import qchem.BasisSet.Atom.Slater.Rk;

void Slater_BS::Register(IBS_Evaluator * eval)
{
    assert(eval);
    eval->Register(&grouper);
}

Rk* Slater_BS::Create(size_t ia,size_t ic,size_t ib,size_t id) const
{
    return new Slater::RkEngine(
        grouper.unique_esv[ia]+grouper.unique_esv[ib],
        grouper.unique_esv[ic]+grouper.unique_esv[id],
        grouper.LMax(ia,ib,ic,id));
}

RVec Slater_BS::loop_4_direct(size_t id, size_t la, size_t lc)  const
{
    const Cacheable* c=Cache4::loop_4(id);
    const Slater::RkEngine* cd = dynamic_cast<const Slater::RkEngine*>(c);
    return cd->Coulomb_Rk(la,lc);
}
RVec Slater_BS::loop_4_exchange(size_t id, size_t la, size_t lc)  const
{
    const Cacheable* c=Cache4::loop_4(id);
    const Slater::RkEngine* cd = dynamic_cast<const Slater::RkEngine*>(c);
    return cd->ExchangeRk(la,lc);
}
    