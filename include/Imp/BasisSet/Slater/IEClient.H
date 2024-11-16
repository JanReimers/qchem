// File: SlaterIEClient
#ifndef _SlaterIEClient_H_
#define _SlaterIEClient_H_

#include <IEClient.H> // <> prevents file search on this Slater folder.
#include "Imp/BasisSet/AtomIEClient.H"

#include "Imp/Integrals/SlaterCD.H"
#include "Imp/BasisSet/BFGrouper.H"
#include "oml/vector.h"

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

//Use the same struct for a full basis set.;
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

    virtual Vector<double> loop_4_direct(size_t id, size_t la, size_t lc)  const
    {
        const Cacheable* c=Cache4::loop_4(es_indices[id-1]);
        const SlaterCD* cd = dynamic_cast<const SlaterCD*>(c);
        return cd->Coulomb_Rk(la,lc);
    }
    virtual Vector<double> loop_4_exchange(size_t id, size_t la, size_t lc)  const
    {
        const Cacheable* c=Cache4::loop_4(es_indices[id-1]);
        const SlaterCD* cd = dynamic_cast<const SlaterCD*>(c);
        return cd->ExchangeRk(la,lc);
    }
}; 

} //namespace Slater

#endif // _SlaterIEClient_H_
