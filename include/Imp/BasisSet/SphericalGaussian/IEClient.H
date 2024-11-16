// File: SphericalGaussianIEClient
#ifndef _SphericalGaussianIEClient_H_
#define _SphericalGaussianIEClient_H_

#include <IEClient.H> // <> prevents file search on this SphericalGaussian folder.
#include "Imp/Integrals/SphericalGaussianCD.H"
#include "Imp/BasisSet/AtomIEClient.H"
#include "Imp/BasisSet/BFGrouper.H"
#include "oml/vector.h"

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
    virtual const SphericalGaussianCD* loop_4(size_t id)  const
    {
        const Cacheable* c=Cache4::loop_4(es_indices[id-1]);
        return dynamic_cast<const SphericalGaussianCD*>(c);
    }
}; 

} //namespace

#endif // _SphericalGaussianIEClient_H_
