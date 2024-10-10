// File SphericalGaussian/BasisSet.H

#include "Imp/BasisSet/SphericalGaussian/BasisSet.H"
#include "Imp/BasisSet/SphericalGaussian/IrrepBasisSet.H"

namespace SphericalGaussian
{


BasisSet::BasisSet(const LinearAlgebraParams& lap,size_t N, double minexp, double maxexp, size_t Lmax)
: BasisGroup() // this makes a integral DB
{
    for (size_t L=0;L<=Lmax;L++)
        Insert(new IrrepBasisSet(lap,GetDataBase(),N,minexp,maxexp,L));
        
}

} //namespace
