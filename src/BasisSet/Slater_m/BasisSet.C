// File Slater_m/BasisSet.H

#include "Imp/BasisSet/Slater_m/BasisSet.H"
#include "Imp/BasisSet/Slater_m/IrrepBasisSet.H"
#include "Imp/BasisSet/Slater_m/IntegralEngine.H"
#include "Imp/BasisSet/SlaterScaler.H"
#include "Imp/Integrals/AngularIntegrals.H"

namespace Slater_m
{


BasisSet::BasisSet(const LAParams& lap,size_t N, double emin, double emax, size_t LMax)
: Slater::BasisSet(new IntegralEngine) // this makes a integral DB
{
    SlaterScaler ss(N,emin,emax,LMax);
    for (size_t L=0;L<=LMax;L++)
        for (int m=-L;m<=(int)L;m++)
            Insert(new IrrepBasisSet(lap,GetDataBase(),this,ss.Get_es(L),L,m));            

}

BasisSet::RVec BasisSet::Coulomb_AngularIntegrals(size_t la, size_t lc, int ma, int mc) const
{
    return AngularIntegrals::Coulomb(la,lc,ma,mc);
}

BasisSet::RVec BasisSet::ExchangeAngularIntegrals(size_t la, size_t lb, int ma, int mb) const
{
    return AngularIntegrals::Exchange(la,lb,ma,mb);
}

} //namespace
