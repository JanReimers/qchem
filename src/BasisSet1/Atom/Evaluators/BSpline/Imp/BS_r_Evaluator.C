// File: BasisSet/Atom/radial/Imp/BSpline_BS_Evaluator.C
module;
#include <cassert>
module BasisSet1.Atom.Evaluators.BSpline.BS;
import qchem.BasisSet1.Atom.BSpline.Rk;


template <size_t K> BSpline_r_BS_Evaluator<K>::~BSpline_r_BS_Evaluator()
{
    delete itsRkCache;
}

template <size_t K> void BSpline_r_BS_Evaluator<K>::Register(IBS_Evaluator * eval)
{
    assert(eval);
    eval->Register(&grouper);
}

template <size_t K> const GLCache* BSpline_r_BS_Evaluator<K>::GetGL(size_t l) const
{
    auto i=grouper.itsGLs.find(l);
    assert(i!=grouper.itsGLs.end());
    return i->second;
}

template <size_t K> void BSpline_r_BS_Evaluator<K>::BuildCache(size_t lmax)
{
    itsRkCache=new BSpline::RkCache_r<K>(grouper.unique_spv,*this->GetGL(lmax),lmax);
}


template <size_t K> Rk* BSpline_r_BS_Evaluator<K>::Create(size_t ia,size_t ic,size_t ib,size_t id) const
{
    assert(itsRkCache);
    // std::cout << "ia,ib,ic,id=" << ia << " " << ib << " " << ic << " " << id << std::endl;
    size_t lmax=grouper.LMax(ia,ib,ic,id);
    const GLCache* gl=this->GetGL(lmax);
    return new BSpline::RkEngine_r(grouper.unique_spv,ia,ib,ic,id,lmax,*gl,*itsRkCache);
}

template <size_t K> double BSpline_r_BS_Evaluator<K>::loop_4_direct(size_t id, size_t la, size_t lc, const rvec11_t& Ak)  const
{
    const Cacheable* c=Cache4::loop_4(id);
    assert(c);
    const BSpline::RkEngine<K>* cd = dynamic_cast<const BSpline::RkEngine<K>*>(c);
    assert(cd);
    return cd->Coulomb_Rk(la,lc,Ak);
}
template <size_t K> double BSpline_r_BS_Evaluator<K>::loop_4_exchange(size_t id, size_t la, size_t lc, const rvec11_t& Ak)  const
{
    const Cacheable* c=Cache4::loop_4(id);
    const BSpline::RkEngine<K>* cd = dynamic_cast<const BSpline::RkEngine<K>*>(c);
    return cd->ExchangeRk(la,lc,Ak);
}

#define INSTANCEk(k) template class BSpline_r_BS_Evaluator<k>;
#include "../Instance.hpp"