// File: SphericalGaussianIEClient.H
#ifndef _SphericalGaussianIEClient_H_
#define _SphericalGaussianIEClient_H_

#include "Imp/BasisSet/AtomIEClient.H"

namespace SphericalGaussian
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
    virtual const Cacheable* Create(size_t ia,size_t ic,size_t ib,size_t id) const;
    virtual Vector<double> loop_4_direct  (size_t id, size_t la, size_t lc)  const;
    virtual Vector<double> loop_4_exchange(size_t id, size_t la, size_t lc)  const;
}; 

} //namespace

#endif // _SphericalGaussianIEClient_H_
