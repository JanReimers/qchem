// File: SlaterIEClient
#ifndef _SlaterIEClient_H_
#define _SlaterIEClient_H_

#include "Imp/BasisSet/AtomIEClient.H"

namespace Slater
{
    
struct IrrepIEClient
    : public virtual ::IrrepIEClient
    , public AtomIrrepIEClient
{
    IrrepIEClient() {};
    IrrepIEClient(size_t N) : AtomIrrepIEClient(N) {};
private:
    virtual double Norm(double e, size_t l) const;
};

struct IEClient 
    : public virtual ::IEClient 
    , public AtomIEClient
{
    
};
} //namespace Slater

#endif // _SlaterIEClient_H_
