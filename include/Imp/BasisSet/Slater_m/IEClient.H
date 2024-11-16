// File: Slater_m/IEClient.H
#ifndef _Slater_m_IEClient_H_
#define _Slater_m_IEClient_H_

#include "Imp/BasisSet/Slater/IEClient.H" 

    
namespace Slater_m
{

typedef Slater::IrrepIEClient IrrepIEClient;    

struct IEClient 
    : public virtual ::IEClient 
    , public AtomIEClient
{
    virtual const Cacheable* Create(size_t ia,size_t ic,size_t ib,size_t id) const;
    virtual Vector<double> loop_4_direct  (size_t id, size_t la, size_t lc)  const;
    virtual Vector<double> loop_4_exchange(size_t id, size_t la, size_t lc)  const;
}; 

} //namespace Slater_m

#endif // _Slater_m_IEClient_H_
