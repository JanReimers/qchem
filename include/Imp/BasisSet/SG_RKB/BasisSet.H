// File: SphericalGaussian_RKB/BasisSet.H
#ifndef _SphericalGaussian_RKB_BasisSet_H_
#define _SphericalGaussian_RKB_BasisSet_H_

#include "IEClient.H"
#include "../TCommon.H"

namespace SphericalGaussian_RKB
{

class DiracBasisSet 
    : public ::BasisSetImp
{
public:
    DiracBasisSet(const LAParams&,size_t N, double minexp, double maxexp, size_t lmax);
    
};

} //namespace

#endif // _SphericalGaussian_RKB_BasisSet_H_
