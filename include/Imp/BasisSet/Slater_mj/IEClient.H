// File: Slater_mj/IEClient.H
#ifndef _Slater_mj_IEClient_H_
#define _Slater_mj_IEClient_H_

#include "Imp/BasisSet/Slater/IEClient.H" 

    
namespace Slater_mj
{
    
typedef Slater::IrrepIEClient Large_IrrepIEClient;  
typedef Slater::IrrepIEClient Small_IrrepIEClient;  

struct Dirac_IrrepIEClient
    : public virtual ::IrrepIEClient
{
    Dirac_IrrepIEClient() : itsLargeIEC(0), itsSmallIEC(0) {};
    virtual size_t size() const;

    void Init(const Large_IrrepIEClient*,const Small_IrrepIEClient*);
    
    const Large_IrrepIEClient* itsLargeIEC;
    const Small_IrrepIEClient* itsSmallIEC;
        
};


} //namespace Slater_mj

#endif // _Slater_mj_IEClient_H_
