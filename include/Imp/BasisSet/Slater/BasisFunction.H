// File: SlaterBF.H  r^l exp(-ar^2) type basis function for an Spherical.
#ifndef _SlaterBF_H_
#define _SlaterBF_H_

#include <BasisFunction.H>

namespace Slater
{
    
class BasisFunction
    : public virtual ::TBasisFunction<double>
{
public:
    BasisFunction();
    BasisFunction(double ex, int n, int l, double norm);
    
    virtual bool   operator==(const ::BasisFunction&) const;

    virtual std::ostream&       Write(std::ostream&) const;
    virtual std::istream&       Read (std::istream&)      ;
    virtual BasisFunction* Clone(        ) const;

    virtual double operator()(const Vec3&) const;
    virtual Vec3   Gradient  (const Vec3&) const;

private:
    double itsExponent;
    int    itsN;
    int    itsL;
    double itsNormalization;
};


} //namespace Slater

#endif //_SlaterBF_H_
