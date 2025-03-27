// File: Atom/ml/Gaussian_BF.H  r^l exp(-ar^2)*Y_lm type basis function.
#ifndef _Atom_ml_Gaussian_BF_H_
#define _Atom_ml_Gaussian_BF_H_

#include <BasisFunction.H>

namespace Atom_ml
{
namespace Gaussian
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
} //namespace

#endif //_Atom_ml_Gaussian_BF_H_
