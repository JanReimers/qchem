// File: Atom/radial/Slater/BS_Common.C  l/ml/kappa/mj independent part of BasisSet for atom Slater basis functions.

#include <vector>
#include <iostream>
#include "radial/Slater/Rk.H"
#include "radial/Slater/BS_Common.H"
import qchem.Irrep_BS;
import oml;

namespace Slater
{
    void BS_Common::Insert(bs_t* bs)
    {
        ::BS_Common::Insert(bs);
        auto iec=dynamic_cast<const IrrepIEClient*>(bs);
        assert(iec);
        Append(iec);
    }
    
    const Cacheable* BS_Common::Create(size_t ia,size_t ic,size_t ib,size_t id) const
    {
        return new Slater::RkEngine(unique_esv[ia]+unique_esv[ib],unique_esv[ic]+unique_esv[id],LMax(ia,ib,ic,id));
    }
    
    
    Vector<double> BS_Common::loop_4_direct(size_t id, size_t la, size_t lc)  const
    {
        const Cacheable* c=Cache4::loop_4(id);
        const Slater::RkEngine* cd = dynamic_cast<const Slater::RkEngine*>(c);
        return cd->Coulomb_Rk(la,lc);
    }
    Vector<double> BS_Common::loop_4_exchange(size_t id, size_t la, size_t lc)  const
    {
        const Cacheable* c=Cache4::loop_4(id);
        const Slater::RkEngine* cd = dynamic_cast<const Slater::RkEngine*>(c);
        return cd->ExchangeRk(la,lc);
    }
    
}