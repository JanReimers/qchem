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
        IrrepBasisSet* ibs=new IrrepBasisSet(lap,GetDataBase(),ss.Get_es(L),L);
        Append(ibs);
        itsIE->Append(ibs); //IECleint
        Insert(ibs);
        
    }
        
    
        
}

} //namespace
