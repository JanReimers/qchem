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
            IrrepBasisSet* ibs=0;
            if (L==0)
                ibs=new IrrepBasisSet(lap,GetDataBase(),NL,ss.emin(L),ss.emax(L),L,m);
            else
            {
                const AtomIrrepIEClient* ibs0=(*this)[1];
                ibs=new IrrepBasisSet(lap,GetDataBase(),NL,ss.emin(L),ibs0->es(NL),L,m);             
            }
            Append(ibs);
            Insert(ibs);            
        }
    }
        
    
        
}

} //namespace
