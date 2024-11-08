// File: SphericalGaussian_m/BasisFunction.H  r^l exp(-ar^2)*Y_lm type basis function.
#ifndef _SphericalGaussian_m_BasisFunction_H_
#define _SphericalGaussian_m_BasisFunction_H_

#include <BasisSet.H>

namespace SphericalGaussian_m
{
    
class BasisFunction
    : public virtual ::TBasisFunction<double>
{
public:
    BasisFunction(                             );
    BasisFunction(double theExponent,int n, int l, int m, double norm);
    
    virtual bool   operator==(const ::BasisFunction&) const;

    virtual std::ostream&       Write(std::ostream&) const;
    virtual std::istream&       Read (std::istream&)      ;
    virtual BasisFunction* Clone(        ) const;

    virtual double operator()(const Vec3&) const;
    virtual Vec3   Gradient  (const Vec3&) const;

private:
    double itsExponent;
    int    n,l,m;
    double itsNormalization;
};


} //namespace

#endif //_SphericalGaussian_m_BasisFunction_H_
