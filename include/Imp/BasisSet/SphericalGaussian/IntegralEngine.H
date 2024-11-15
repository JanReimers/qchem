// File: SphericalGaussianIE1.H  Integral Engine for Spherical gaussians.
#ifndef _SphericalGaussianIE1_H_
#define _SphericalGaussianIE1_H_


#include <AnalyticIE.H>
#include "Imp/BasisSet/AtomIE.H"
#include "oml/vector.h"


namespace SphericalGaussian
{
    
class IrrepIEClient;
    
class IntegralEngine
    : public virtual AnalyticIE<double>
    , public AtomIE
{
public:
    virtual void Make4C(ERI4& J, ERI4& K, const ::IEClient*) const; //Used for Hartree-Fock


private:
    virtual double Overlap(double ea, double eb,size_t l) const;
    virtual double Kinetic(double ea ,double eb,size_t l) const;
    virtual double Nuclear(double ea, double eb,size_t l) const;
    virtual double Charge (double ea,           size_t l) const;
    virtual double Repulsion(double eab, double ec,size_t la,size_t lc) const;
    
    static double FourPi2;
};

} //namespace

#endif //_SphericalGaussianIE1_H_
