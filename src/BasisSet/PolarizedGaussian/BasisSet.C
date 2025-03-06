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
    Insert(new IrrepBasisSet(lap,GetDataBase(),this,reader,cl));
}

BasisSet::BasisSet(const LAParams& lap, size_t N, double emin, double emax, size_t LMax, const Cluster* cl)
: BasisSetImp(new IntegralEngine) // this makes a integral DB
{
    GaussianScaler gs(N,emin,emax,LMax);
    if (cl->GetNumAtoms()>1)
        Insert(new IrrepBasisSet(lap,GetDataBase(),this,gs.Get_es(0),LMax,cl));        
    else
    {
        assert(cl->GetNumAtoms()==1);
        for (size_t L=0;L<=LMax;L++)
            Insert(new IrrepBasisSet(lap,GetDataBase(),this,gs.Get_es(L),L));  
            
    }
}

void BasisSet::Insert(bs_t* bs)
{
    BasisSetImp::Insert(bs);
    Append(bs);
}
} //namespace
