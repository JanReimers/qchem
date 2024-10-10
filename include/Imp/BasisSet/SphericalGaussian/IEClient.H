// File: SphericalGaussianIEClient
#ifndef _SphericalGaussianIEClient_H_
#define _SphericalGaussianIEClient_H_

#include <IEClient.H> // <> prevents file search on this SphericalGaussian folder.
#include "oml/vector.h"

namespace SphericalGaussian
{
    
struct IEData
{
    IEData() {}; 
    IEData(size_t N) : Ls(N), es(N), ns(N) {};
    Vector<size_t> Ls;  //Angular momentum.
    Vector<double> es; //The orbital exponents.
    Vector<double> ns; //Normalization constants
};

    
struct IrrepIEClient
    : public virtual ::IrrepIEClient
    , public IEData
{
    IrrepIEClient() {};
    IrrepIEClient(size_t N) : IEData(N) {};
    void Init(double minexp,double maxexp,size_t L);
    
    virtual size_t size() const {return es.size();}
    typedef std::tuple<int,double,double> bf_tuple;
    bf_tuple operator()(index_t i) const {return std::make_tuple(Ls(i),es(i),ns(i));}
};

//Use the same struct for a full basis set.;
struct IEClient 
    : public virtual ::IEClient 
    , public IEData
{
    IEClient() {};
    IEClient(size_t N) : IEData(N) {};
    virtual size_t size() const {return es.size();}
}; 

} //namespace

#endif // _SphericalGaussianIEClient_H_
