// File: BasisSet/Atom/radial/Imp/BSpline_BS.C
module;
#include <cassert>
module BasisSet.Atom.BSpline_BS;
import qchem.BasisSet.Atom.Internal.radial.BSpline.Rk;
       

template <size_t K> void BSpline_BS<K>::Register(IBS_Evaluator * eval)
{
    assert(eval);
    eval->Register(&grouper);
}

template <size_t K> Rk* BSpline_BS<K>::Create(size_t ia,size_t ic,size_t ib,size_t id) const
{
    assert(false);
    return 0;
    // return new BSpline::RkEngine(
    //     grouper.unique_esv[ia]+grouper.unique_esv[ib],
    //     grouper.unique_esv[ic]+grouper.unique_esv[id],
    //     grouper.LMax(ia,ib,ic,id));
}

template <size_t K> RVec BSpline_BS<K>::loop_4_direct(size_t id, size_t la, size_t lc)  const
{
    const Cacheable* c=Cache4::loop_4(id);
    const BSpline::RkEngine<6>* cd = dynamic_cast<const BSpline::RkEngine<6>*>(c);
    return cd->Coulomb_Rk(la,lc);
}
template <size_t K> RVec BSpline_BS<K>::loop_4_exchange(size_t id, size_t la, size_t lc)  const
{
    const Cacheable* c=Cache4::loop_4(id);
    const BSpline::RkEngine<6>* cd = dynamic_cast<const BSpline::RkEngine<6>*>(c);
    return cd->ExchangeRk(la,lc);
}

#define INSTANCEk(k) template class BSpline_BS<k>;
#include "../BSpline/Instance.hpp"