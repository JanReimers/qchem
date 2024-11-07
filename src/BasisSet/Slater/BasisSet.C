// File Slater/BasisSet.H

#include "Imp/BasisSet/Slater/BasisSet.H"
#include "Imp/BasisSet/Slater/IrrepBasisSet.H"
#include "Imp/BasisSet/Slater/IntegralEngine.H"
#include "Imp/BasisSet/SlaterScaler.H"

namespace Slater
{


BasisSet::BasisSet(const LAParams& lap,size_t N, double emin, double emax, size_t LMax)
: BasisSetImp(new IntegralEngine) // this makes a integral DB
{
    SlaterScaler ss(N,emin,emax,LMax);
    for (size_t L=0;L<=LMax;L++)
    {
        size_t  NL=ss.N(L);
        IrrepBasisSet* ibs= 
            L==0 ?
                new IrrepBasisSet(lap,GetDataBase(),NL,ss.emin(L),ss.emax(L),L)
            :
                new IrrepBasisSet(lap,GetDataBase(),NL,ss.emin(L),es(NL),L); 
        Append(ibs);
        Insert(ibs);
        
    }
        
    
        
}

} //namespace
