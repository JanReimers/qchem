// File: SphericalGaussian_m/IEClient.H
#ifndef _SphericalGaussian_m_IEClient_H_
#define _SphericalGaussian_m_IEClient_H_

#include "Imp/BasisSet/SphericalGaussian/IEClient.H" 
#include "Imp/Integrals/SphericalGaussianCD.H"
#include "Imp/BasisSet/BFGrouper.H"
#include "oml/vector.h"

namespace SphericalGaussian_m
{
    
typedef SphericalGaussian::IrrepIEClient IrrepIEClient;

struct IEClient 
    : public virtual ::IEClient 
    , public AtomIEClient
{
    virtual const Cacheable* Create(size_t ia,size_t ic,size_t ib,size_t id) const;
    
    virtual const SphericalGaussianCD* loop_4(size_t id)  const
    {
        const Cacheable* c=Cache4::loop_4(es_indices[id-1]);
        return dynamic_cast<const SphericalGaussianCD*>(c);
    }



    
}; 

} //namespace SphericalGaussian_m

#endif // _Slater_m_IEClient_H_
