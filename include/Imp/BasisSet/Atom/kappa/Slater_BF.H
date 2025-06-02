// File: Atom/kappa/Slater_BF.H  Slater basis function with Restricted Kinetic Balance (RKB).
#ifndef _Atom_kappa_Slater_BF_H_
#define _Atom_kappa_Slater_BF_H_

#include <BasisFunction.H>

namespace Atom_kappa
{
namespace Slater
{

class Small_BasisFunction;
    
class Large_BasisFunction
    : public virtual ::TBasisFunction<double>
{
public:
    Large_BasisFunction() {};
    Large_BasisFunction(double ex, int kappa, int mj, double norm);
    
    virtual std::ostream&    Write(std::ostream&) const;
    virtual ::BasisFunction* Clone(        ) const;

    virtual double operator()(const Vec3&) const;
    virtual Vec3   Gradient  (const Vec3&) const;
private:
    friend class Small_BasisFunction;
    double itsExponent;
    int kappa,mj,l;
    double itsNormalization;
};

//
//  Derived from the LargeBF P(r)=r^l*exp(-e*r) as:
//
//                                  -e*r^l*exp(-e*r), kappa<0
//     Q(r)=(d/dr+(1+kappa)/r)g = {
//                                  ((2l+1)/r-e)*r^lexp(-e*r), kappa>0
//
class Small_BasisFunction
    : public virtual ::TBasisFunction<double>
{
public:
    Small_BasisFunction();
    Small_BasisFunction(const Large_BasisFunction*,double norm);
    
    virtual std::ostream&    Write(std::ostream&) const;
    virtual ::BasisFunction* Clone(        ) const;

    virtual double operator()(const Vec3&) const;
    virtual Vec3   Gradient  (const Vec3&) const;
private:
    const Large_BasisFunction* Pr;
    double itsNormalization;
};

}} //namespace _Atom_kappa::Slater_

#endif //_Atom_kappa_Slater_BF_H_
