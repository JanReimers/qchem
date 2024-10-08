// File: SphericalGaussianIEClient
#ifndef _SphericalGaussianIEClient_H_
#define _SphericalGaussianIEClient_H_

#include <IEClient.H> // <> prevents file search on this SphericalGaussian folder.
#include "oml/vector.h"

struct SphericalGaussianIEClient
    : public virtual IrrepIEClient
{
    SphericalGaussianIEClient() {};
    SphericalGaussianIEClient(size_t N) : Ls(N), es(N), ns(N) {};
    void Init(double minexp,double maxexp,size_t L);
    
    virtual size_t size() const {return es.size();}

    Vector<size_t> Ls;  //Angular momentum.
    Vector<double> es; //The orbital exponents.
    Vector<double> ns; //Normalization constants
};
   

#endif // _SphericalGaussianIEClient_H_
