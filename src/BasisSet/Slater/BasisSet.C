// File Slater/BasisSet.H

#include "Imp/BasisSet/Slater/BasisSet.H"
#include "Imp/BasisSet/Slater/IrrepBasisSet.H"
#include "Imp/BasisSet/Slater/IntegralEngine.H"
#include "Imp/BasisSet/SlaterScaler.H"
#include "Imp/Integrals/SlaterCD.H"
#include "Imp/Integrals/SlaterIntegrals.H"
#include "Imp/Integrals/AngularIntegrals.H"

namespace Slater
{


BasisSet::BasisSet(const LAParams& lap,size_t N, double emin, double emax, size_t LMax)
: BasisSetImp(new IntegralEngine) // this makes a integral DB
{
    SlaterScaler ss(N,emin,emax,LMax);
    for (size_t L=0;L<=LMax;L++)
        Insert(new Orbital_IBS(lap,this,ss.Get_es(L),L));
        
}

void BasisSet::Insert(bs_t* bs)
{
    BasisSetImp::Insert(bs);
    Append(bs);
}

BasisSet::RVec BasisSet::Coulomb_AngularIntegrals(size_t la, size_t lc, int, int) const
{
    return AngularIntegrals::Coulomb(la,lc);
}

BasisSet::RVec BasisSet::ExchangeAngularIntegrals(size_t la, size_t lb, int, int) const
{
    return AngularIntegrals::Exchange(la,lb);
}

const Cacheable* BasisSet::Create(size_t ia,size_t ic,size_t ib,size_t id) const
{
//        cout << "new " << ia << " " << ib << " " << ic << " " << id << endl;
//        cout << "new " << unique_esv[ia] << " " << unique_esv[ib] << " " << unique_esv[ic] << " " << unique_esv[id] << endl;
    return new SlaterCD(unique_esv[ia]+unique_esv[ib],unique_esv[ic]+unique_esv[id],LMax(ia,ib,ic,id));
}


Vector<double> BasisSet::loop_4_direct(size_t id, size_t la, size_t lc)  const
{
    const Cacheable* c=Cache4::loop_4(id);
    const SlaterCD* cd = dynamic_cast<const SlaterCD*>(c);
    return cd->Coulomb_Rk(la,lc);
}
Vector<double> BasisSet::loop_4_exchange(size_t id, size_t la, size_t lc)  const
{
    const Cacheable* c=Cache4::loop_4(id);
    const SlaterCD* cd = dynamic_cast<const SlaterCD*>(c);
    return cd->ExchangeRk(la,lc);
}


} //namespace
