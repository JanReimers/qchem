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

    virtual Vector<double> loop_4_direct(size_t id, size_t la, size_t lc)  const
    {
        const Cacheable* c=Cache4::loop_4(es_indices[id-1]);
        const SphericalGaussianCD* cd = dynamic_cast<const SphericalGaussianCD*>(c);
        return cd->Coulomb_Rk(la,lc);
    }
    virtual Vector<double> loop_4_exchange(size_t id, size_t la, size_t lc)  const
    {
        const Cacheable* c=Cache4::loop_4(es_indices[id-1]);
        const SphericalGaussianCD* cd = dynamic_cast<const SphericalGaussianCD*>(c);
        return cd->ExchangeRk(la,lc);
    }

    
}; 

} //namespace SphericalGaussian_m

#endif // _Slater_m_IEClient_H_
