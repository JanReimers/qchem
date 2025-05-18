// File: Atom/radial/BSpline/BS_Common.H  l/ml/kappa/mj independent part of BasisSet for atom BSpline Basis Sets.

#include "Imp/BasisSet/Atom/radial/BSpline/BS_Common.H"
#include "Imp/BasisSet/Atom/radial/BSpline/Rk.H"
#include <Irrep_BS.H>
#include "oml/vector.h"

namespace BSpline
{
    template <size_t K> void BS_Common<K>::Insert(bs_t* bs)
    {
        ::BS_Common::Insert(bs);
        this->Append(bs); 
    }
    
    template <size_t K> const Cacheable* BS_Common<K>::Create(size_t ia,size_t ic,size_t ib,size_t id) const
    {
        // std::cout << "ia,ib,ic,id=" << ia << " " << ib << " " << ic << " " << id << std::endl;
        size_t lmax=LMax(ia,ib,ic,id);
        const GLCache* gl=this->GetGL(lmax);
        return new BSpline::RkEngine(unique_spv[ia],unique_spv[ib],unique_spv[ic],unique_spv[id],lmax,*gl);

    }
    
    
    template <size_t K> Vector<double> BS_Common<K>::loop_4_direct(size_t id, size_t la, size_t lc)  const
    {
        const Cacheable* c=Cache4::loop_4(id);
        const BSpline::RkEngine<K>* cd = dynamic_cast<const BSpline::RkEngine<K>*>(c);
        return cd->Coulomb_Rk(la,lc);
    }
    template <size_t K> Vector<double> BS_Common<K>::loop_4_exchange(size_t id, size_t la, size_t lc)  const
    {
        const Cacheable* c=Cache4::loop_4(id);
        const BSpline::RkEngine<K>* cd = dynamic_cast<const BSpline::RkEngine<K>*>(c);
        return cd->ExchangeRk(la,lc);
    }
    
    template class BS_Common<6>;
}