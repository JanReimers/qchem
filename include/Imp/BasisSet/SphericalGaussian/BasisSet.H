// File SphericalGaussian/BasisSet.H
#ifndef _SphericalGaussian_BasisSet_H_
#define _SphericalGaussian_BasisSet_H_

#include <BasisSet.H>

namespace SphericalGaussian
{

class BasisSet : public BasisGroup
{
public:
    BasisSet(const LinearAlgebraParams&,size_t N, double minexp, double maxexp, size_t Lmax);
    
};

} //namespace

#endif // _SphericalGaussian_BasisSet_H_
