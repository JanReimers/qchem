// File Slater_m/BasisSet.H

#include "Imp/BasisSet/Slater_m/BasisSet.H"
#include "Imp/BasisSet/Slater_m/IrrepBasisSet.H"
#include "Imp/BasisSet/Slater_m/IntegralEngine.H"
#include "Imp/BasisSet/SlaterScaler.H"

namespace Slater_m
{


BasisSet::BasisSet(const LAParams& lap,size_t N, double emin, double emax, size_t LMax)
: BasisSetImp(new IntegralEngine) // this makes a integral DB
{
    SlaterScaler ss(N,emin,emax,LMax);
    for (size_t L=0;L<=LMax;L++)
        for (int m=-L;m<=(int)L;m++)
            Insert(new IrrepBasisSet(lap,GetDataBase(),ss.Get_es(L),L,m));            

}

} //namespace
