// File PolarizedGaussian/BasisSet.H
#ifndef _PolarizedGaussian_BasisSet_H_
#define _PolarizedGaussian_BasisSet_H_

#include "IEClient.H"
#include "../TCommon.H"

namespace PolarizedGaussian
{

class Reader;

class BasisSet 
    : public virtual IEClient
    , public BasisSetImp
{
public:
    BasisSet() {};
    BasisSet(const LAParams&, Reader*, const Cluster*);
    //! Make a primatives dentered at R=(0.0,0) with polarizations up to LMax
    BasisSet(const LAParams&, size_t N, double emin, double emax, size_t LMax);
    
};

} //namespace

#endif // _PolarizedGaussian

