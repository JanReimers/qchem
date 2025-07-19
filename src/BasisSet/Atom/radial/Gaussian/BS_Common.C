// File: Atom/radial/Gaussian/BS_Common.C  l/ml/kappa/mj independent part of BasisSet for atom Gaussians.

#include <vector>
#include <iostream>
#include <cassert>
#include "radial/Gaussian/Rk.H"
#include "radial/Gaussian/BS_Common.H"
import qchem.Irrep_BS;
import oml;

namespace Gaussian
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
    return new ::Gaussian::RkEngine(unique_esv[ia]+unique_esv[ib],unique_esv[ic]+unique_esv[id],LMax(ia,ib,ic,id));
}

Vector<double>  BS_Common::loop_4_direct(size_t id, size_t la, size_t lc)  const
{
    const Cacheable* c=Cache4::loop_4(id);
    const ::Gaussian::RkEngine* cd = dynamic_cast<const ::Gaussian::RkEngine*>(c);
    return cd->Coulomb_Rk(la,lc);
}
Vector<double>  BS_Common::loop_4_exchange(size_t id, size_t la, size_t lc)  const
{
    const Cacheable* c=Cache4::loop_4(id);
    const ::Gaussian::RkEngine* cd = dynamic_cast<const ::Gaussian::RkEngine*>(c);
    return cd->ExchangeRk(la,lc);
}

} //namespace