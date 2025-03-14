// File Slater/BasisSet.H

#include "Imp/BasisSet/Slater/BasisSet.H"
#include "Imp/BasisSet/Slater/IrrepBasisSet.H"
#include "Imp/BasisSet/SlaterScaler.H"
#include "Imp/Integrals/SlaterCD.H"

namespace Slater
{


BasisSet::BasisSet(const LAParams& lap,size_t N, double emin, double emax, size_t LMax)
{
    SlaterScaler ss(N,emin,emax,LMax);
    for (size_t L=0;L<=LMax;L++)
        Insert(new Orbital_IBS(lap,this,ss.Get_es(L),L));
        
}

void BS_Common::Insert(bs_t* bs)
{
    BasisSetImp::Insert(bs);
    Append(bs);
}

const Cacheable* BS_Common::Create(size_t ia,size_t ic,size_t ib,size_t id) const
{
//        cout << "new " << ia << " " << ib << " " << ic << " " << id << endl;
//        cout << "new " << unique_esv[ia] << " " << unique_esv[ib] << " " << unique_esv[ic] << " " << unique_esv[id] << endl;
    return new SlaterCD(unique_esv[ia]+unique_esv[ib],unique_esv[ic]+unique_esv[id],LMax(ia,ib,ic,id));
}


Vector<double> BS_Common::loop_4_direct(size_t id, size_t la, size_t lc)  const
{
    const Cacheable* c=Cache4::loop_4(id);
    const SlaterCD* cd = dynamic_cast<const SlaterCD*>(c);
    return cd->Coulomb_Rk(la,lc);
}
Vector<double> BS_Common::loop_4_exchange(size_t id, size_t la, size_t lc)  const
{
    const Cacheable* c=Cache4::loop_4(id);
    const SlaterCD* cd = dynamic_cast<const SlaterCD*>(c);
    return cd->ExchangeRk(la,lc);
}


} //namespace
