// File: Guassian/BasisFunction.H  r^l exp(-a*r^2) type Gaussian basis function.
#ifndef _GaussianBF_H_
#define _GaussianBF_H_

#include <BasisFunction.H>

namespace Gaussian
{
    
class BasisFunction
    : public virtual ::TBasisFunction<double>
{
public:
    BasisFunction(                             );
    BasisFunction(double theExponent, int theL, double norm);
    
    virtual bool   operator==(const ::BasisFunction&) const;

    virtual std::ostream&       Write(std::ostream&) const;
    virtual std::istream&       Read (std::istream&)      ;
    virtual BasisFunction* Clone(        ) const;

    virtual double operator()(const Vec3&) const;
    virtual Vec3   Gradient  (const Vec3&) const;

private:
    double itsExponent;
    int    itsL;
    double itsNormalization;
};


} //namespace

#endif //_GaussianBF_H_
