// File PolarizedGaussian/BasisSet.H
#ifndef _PolarizedGaussian_BasisSet_H_
#define _PolarizedGaussian_BasisSet_H_

#include "IEClient.H"
#include <BasisSet.H>

namespace PolarizedGaussian
{

class Reader;

class BasisSet 
    : public virtual IEClient
    , public BasisGroup
{
public:
    BasisSet() {};
    BasisSet(const LinearAlgebraParams&, Reader*, const Cluster*);
    
};

} //namespace

#endif // _PolarizedGaussian

