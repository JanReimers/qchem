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

struct IEData
{
    IEData() {}; 
    IEData(size_t N) : Ls(N), es(N), ns(N) {};
    Vector<size_t> Ls;  //Angular momentum.
    Vector<double> es; //The orbital exponents.
    Vector<double> ns; //Normalization constants
};


//Use the same struct for a full basis set.;
struct IEClient 
    : public virtual ::IEClient 
    , public virtual Cache4
    , public IEData
    , public BFGrouper
{
    IEClient() {};
    IEClient(size_t N) : IEData(N) {};
    void Append(const IrrepIEClient*);
    
    virtual size_t size() const {return es.size();}
    size_t LMax() const {return Max(Ls);}
    virtual const Cacheable* Create(size_t ia,size_t ic,size_t ib,size_t id) const;
    void loop_1(size_t ia) const {Cache4::loop_1(es_indices[ia-1]);}
    void loop_2(size_t ic) const {Cache4::loop_2(es_indices[ic-1]);}
    void loop_3(size_t ib) const {Cache4::loop_3(es_indices[ib-1]);}

    virtual const SphericalGaussianCD* loop_4(size_t id)  const
    {
        const Cacheable* c=Cache4::loop_4(es_indices[id-1]);
        return dynamic_cast<const SphericalGaussianCD*>(c);
    }
}; 

} //namespace

#endif // _SphericalGaussianIEClient_H_
