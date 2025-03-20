// File: SphericalGaussian_m/BasisSet.H
#ifndef _SphericalGaussian_m_BasisSet_H_
#define _SphericalGaussian_m_BasisSet_H_

#include "IEClient.H"
#include "../TCommon.H"
#include "Imp/BasisSet/SphericalGaussian/BasisSet.H"

namespace PolarizedGaussian
{
    class Reader;    
}

namespace SphericalGaussian_m
{

class BasisSet 
: public SphericalGaussian::BS_Common
, public ::AtomIE_BS_2E_lm //Pick angular integrals.
{
public:
    BasisSet() {};
    BasisSet(const LAParams&,size_t N, double minexp, double maxexp, size_t Lmax);
    BasisSet(const LAParams&,PolarizedGaussian::Reader*, const Atom*);
};

} //namespace SphericalGaussian_m

#endif // _SphericalGaussian_m_BasisSet_H_
