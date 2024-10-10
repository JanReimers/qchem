// File PolarizedGaussian/BasisSet.H

#include "Imp/BasisSet/PolarizedGaussian/BasisSet.H"
#include "Imp/BasisSet/PolarizedGaussian/IrrepBasisSet.H"

namespace PolarizedGaussian
{


BasisSet::BasisSet(const LinearAlgebraParams& lap, Reader* reader, const Cluster* cl)
: BasisGroup() // this makes a integral DB
{
    Insert(new IrrepBasisSet(lap,GetDataBase(),reader,cl));
}

} //namespace
