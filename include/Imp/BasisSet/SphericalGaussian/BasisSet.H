// File SphericalGaussian/BasisSet.H

#include <BasisSet.H>

namespace SphericalGaussian
{

class BasisSet : public BasisGroup
{
public:
    BasisSet(const LinearAlgebraParams&,size_t N, double minexp, double maxexp, size_t Lmax);
    
};

} //namespace
