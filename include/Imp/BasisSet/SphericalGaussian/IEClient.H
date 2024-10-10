// File: SphericalGaussianIEClient
#ifndef _SphericalGaussianIEClient_H_
#define _SphericalGaussianIEClient_H_

#include <IEClient.H> // <> prevents file search on this SphericalGaussian folder.
#include "oml/vector.h"

namespace SphericalGaussian
{
    
struct IrrepIEClient
    : public virtual ::IrrepIEClient
{
    IrrepIEClient() {};
    IrrepIEClient(size_t N) : Ls(N), es(N), ns(N) {};
    void Init(double minexp,double maxexp,size_t L);
    
    virtual size_t size() const {return es.size();}
    typedef std::tuple<int,double,double> bf_tuple;
    bf_tuple operator()(index_t i) const {return std::make_tuple(Ls(i),es(i),ns(i));}

    Vector<size_t> Ls;  //Angular momentum.
    Vector<double> es; //The orbital exponents.
    Vector<double> ns; //Normalization constants
};

} //namespace

#endif // _SphericalGaussianIEClient_H_
