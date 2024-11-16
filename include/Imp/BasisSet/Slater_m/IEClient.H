// File: Slater_m/IEClient.H
#ifndef _Slater_m_IEClient_H_
#define _Slater_m_IEClient_H_

#include "Imp/BasisSet/Slater/IEClient.H" 
#include "Imp/Integrals/SlaterCD.H"
#include "Imp/BasisSet/BFGrouper.H"
    
namespace Slater_m
{

typedef Slater::IrrepIEClient IrrepIEClient;    

struct IEClient 
    : public virtual ::IEClient 
    , public AtomIEClient
{
    virtual const Cacheable* Create(size_t ia,size_t ic,size_t ib,size_t id) const;
    
    virtual const SlaterCD* loop_4(size_t id)  const
    {
        const Cacheable* c=Cache4::loop_4(es_indices[id-1]);
        return dynamic_cast<const SlaterCD*>(c);
    }
 
}; 

} //namespace Slater_m

#endif // _Slater_m_IEClient_H_
