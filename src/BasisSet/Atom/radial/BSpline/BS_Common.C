// File: Atom/radial/BSpline/BS_Common.H  l/ml/kappa/mj independent part of BasisSet for atom BSpline Basis Sets.

#include "Imp/BasisSet/Atom/radial/BSpline/BS_Common.H"
#include "Imp/BasisSet/Atom/radial/BSpline/Rk.H"
#include <Irrep_BS.H>
#include "oml/vector.h"

namespace BSpline
{
    void BS_Common::Insert(bs_t* bs)
    {
        ::BS_Common::Insert(bs);
        //Append(bs); BFGrouper not set up yet
    }
    
    const Cacheable* BS_Common::Create(size_t ia,size_t ic,size_t ib,size_t id) const
    {
        // return new BSpline::RkEngine(unique_esv[ia]+unique_esv[ib],unique_esv[ic]+unique_esv[id],LMax(ia,ib,ic,id));
        return 0;
    }
    
    
    Vector<double> BS_Common::loop_4_direct(size_t id, size_t la, size_t lc)  const
    {
        // const Cacheable* c=Cache4::loop_4(id);
        // const BSpline::RkEngine* cd = dynamic_cast<const BSpline::RkEngine*>(c);
        // return cd->Coulomb_Rk(la,lc);
        return Vector<double>();
    }
    Vector<double> BS_Common::loop_4_exchange(size_t id, size_t la, size_t lc)  const
    {
        // const Cacheable* c=Cache4::loop_4(id);
        // const BSpline::RkEngine* cd = dynamic_cast<const BSpline::RkEngine*>(c);
        // return cd->ExchangeRk(la,lc);
        return Vector<double>();
    }
    
}