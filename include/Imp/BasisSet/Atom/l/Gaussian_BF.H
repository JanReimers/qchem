// File: Atom/l/Gaussian_BF.H r^l exp(-a*r^2) type Gaussian basis function.
#ifndef _Atom_l_Gaussian_BF_H_
#define _Atom_l_Gaussian_BF_H_

#include <BasisSet/BasisFunction.H>

namespace Atoml
{
namespace Gaussian
{
    
class BasisFunction
    : public virtual ::TBasisFunction<double>
{
public:
    BasisFunction(                             );
    BasisFunction(double theExponent, int theL, double norm);
    
    virtual std::ostream&  Write(std::ostream&) const;
    virtual BasisFunction* Clone(        ) const;

    virtual double operator()(const Vec3&) const;
    virtual Vec3   Gradient  (const Vec3&) const;

private:
    double itsExponent;
    int    itsL;
    double itsNormalization;
};

}} //namespace

#endif //_Atom_l_Gaussian_BF_H_
