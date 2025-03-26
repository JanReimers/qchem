// File SphericalGaussian/BasisSet.H

#include "Imp/BasisSet/SphericalGaussian/BasisSet.H"
#include "Imp/BasisSet/SphericalGaussian/IrrepBasisSet.H"
#include "Imp/BasisSet/Atom/radial/Gaussian/ExponentScaler.H"
#include "Imp/BasisSet/Atom/radial/Gaussian/Rk.H"
#include "Imp/Integrals/AngularIntegrals.H"

namespace SphericalGaussian
{


BasisSet::BasisSet(size_t N, double emin, double emax, size_t LMax)
{
    Gaussian::ExponentScaler gs(N,emin,emax,LMax);
    for (size_t L=0;L<=LMax;L++)
        Insert(new Orbital_IBS(this,gs.Get_es(L),L)); //Common with optr_vector     
}

void BS_Common::Insert(bs_t* bs)
{
    BasisSetImp::Insert(bs);
    Append(bs);
}

const Cacheable* BS_Common::Create(size_t ia,size_t ic,size_t ib,size_t id) const
{
    return new Gaussian::RkEngine(unique_esv[ia]+unique_esv[ib],unique_esv[ic]+unique_esv[id],LMax(ia,ib,ic,id));
}

Vector<double>  BS_Common::loop_4_direct(size_t id, size_t la, size_t lc)  const
{
    const Cacheable* c=Cache4::loop_4(id);
    const Gaussian::RkEngine* cd = dynamic_cast<const Gaussian::RkEngine*>(c);
    return cd->Coulomb_Rk(la,lc);
}
Vector<double>  BS_Common::loop_4_exchange(size_t id, size_t la, size_t lc)  const
{
    const Cacheable* c=Cache4::loop_4(id);
    const Gaussian::RkEngine* cd = dynamic_cast<const Gaussian::RkEngine*>(c);
    return cd->ExchangeRk(la,lc);
}


} //namespace
