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
    for (int L=0;L<=LMax;L++)
    {
        size_t  NL=ss.N(L);
        for (int m=-L;m<=L;m++)
        {
            IrrepBasisSet* ibs=new IrrepBasisSet(lap,GetDataBase(),ss.Get_es(L),L,m);
            Append(ibs);
            Insert(ibs);            
        }
    }
        
    
        
}

} //namespace
