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
    
}; 

} //namespace Slater_m

#endif // _Slater_m_IEClient_H_
