// File PolarizedGaussian/BasisSet.H

#include "Imp/BasisSet/PolarizedGaussian/BasisSet.H"
#include "Imp/BasisSet/PolarizedGaussian/IrrepBasisSet.H"
#include "Imp/BasisSet/PolarizedGaussian/IntegralEngine.H"
#include "Imp/BasisSet/GaussianScaler.H"
#include <Cluster.H>

namespace PolarizedGaussian
{


BasisSet::BasisSet(const LAParams& lap, Reader* reader, const Cluster* cl)
: BasisSetImp(new IntegralEngine) // this makes a integral DB
{
    IrrepBasisSet* ibs=new IrrepBasisSet(lap,GetDataBase(),reader,cl);
    Append(ibs);
    Insert(ibs);
}

BasisSet::BasisSet(const LAParams& lap, size_t N, double emin, double emax, size_t LMax, const Cluster* cl)
: BasisSetImp(new IntegralEngine) // this makes a integral DB
{
    GaussianScaler gs(N,emin,emax,LMax);
    if (cl->GetNumAtoms()>1)
    {
        IrrepBasisSet* ibs=new IrrepBasisSet(lap,GetDataBase(),gs.Get_es(0),LMax,cl);
        Append(ibs);
        Insert(ibs);        
    }
    else
    {
        assert(cl->GetNumAtoms()==1);
        for (size_t L=0;L<=LMax;L++)
        {
            IrrepBasisSet* ibs=new IrrepBasisSet(lap,GetDataBase(),gs.Get_es(L),L);
            Append(ibs);
            Insert(ibs);  
        }
            
    }
}
} //namespace
