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

} //namespace

#endif // _SphericalGaussianIEClient_H_
