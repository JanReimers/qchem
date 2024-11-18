// File: SphericalGaussian_m/IEClient.H
#ifndef _SphericalGaussian_m_IEClient_H_
#define _SphericalGaussian_m_IEClient_H_

#include "Imp/BasisSet/SphericalGaussian/IEClient.H" 

namespace SphericalGaussian_m
{
    
typedef SphericalGaussian::IrrepIEClient IrrepIEClient;

struct IEClient 
    : public virtual ::IEClient 
    , public AtomIEClient
{
  
}; 

} //namespace SphericalGaussian_m

#endif // _Slater_m_IEClient_H_
