// File: SphericalGaussian_m/BasisSet.H
#ifndef _SphericalGaussian_m_BasisSet_H_
#define _SphericalGaussian_m_BasisSet_H_

#include "IEClient.H"
#include "../TCommon.H"

namespace SphericalGaussian_m
{

class BasisSet 
    : public SphericalGaussian_m::IEClient
    , public ::BasisSetImp
{
public:
    BasisSet() {};
    BasisSet(const LAParams&,size_t N, double minexp, double maxexp, size_t Lmax);
    
};

} //namespace SphericalGaussian_m

#endif // _SphericalGaussian_m_BasisSet_H_
