// File Slater/BasisSet.H

#include "Imp/BasisSet/Slater/BasisSet.H"
#include "Imp/BasisSet/Slater/IrrepBasisSet.H"
#include "Imp/BasisSet/Slater/IntegralEngine.H"

namespace Slater
{


BasisSet::BasisSet(const LAParams& lap,size_t N, double minexp, double maxexp, size_t Lmax)
: BasisSetImp(new IntegralEngine) // this makes a integral DB
{
    for (size_t L=0;L<=Lmax;L++)
    {
        IrrepBasisSet* ibs=new IrrepBasisSet(lap,GetDataBase(),N,minexp,maxexp,L);
        Append(ibs);
        Insert(ibs);
        
    }
        
    
        
}

} //namespace
