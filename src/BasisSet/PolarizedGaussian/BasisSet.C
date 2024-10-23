// File PolarizedGaussian/BasisSet.H

#include "Imp/BasisSet/PolarizedGaussian/BasisSet.H"
#include "Imp/BasisSet/PolarizedGaussian/IrrepBasisSet.H"
#include "Imp/BasisSet/PolarizedGaussian/IntegralEngine.H"

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
    IrrepBasisSet* ibs=new IrrepBasisSet(lap,GetDataBase(),N,emin,emax,LMax,cl);
    Append(ibs);
    Insert(ibs);
}
} //namespace
