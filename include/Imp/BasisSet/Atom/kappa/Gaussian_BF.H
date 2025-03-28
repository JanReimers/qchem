// File: Atom/kappa/Gaussian_BF.H  Gaussians with Restricted Kinetic Balance (RKB).
#ifndef _Atom_kappa_Gaussian_BF_H_
#define _Atom_kappa_Gaussian_BF_H_

#include <BasisFunction.H>

namespace Atom_kappa
{
namespace Gaussian
{

class Small_BasisFunction;

class Large_BasisFunction
    : public virtual ::TBasisFunction<double>
{
public:
    Large_BasisFunction(                             );
    Large_BasisFunction(double theExponent, int kappa, double norm);
    
    virtual bool   operator==(const ::BasisFunction&) const;

    virtual std::ostream&       Write(std::ostream&) const;
    virtual std::istream&       Read (std::istream&)      ;
    virtual BasisFunction* Clone(        ) const;

    virtual double operator()(const Vec3&) const;
    virtual Vec3   Gradient  (const Vec3&) const;

private:
    friend class Small_BasisFunction;

    double itsExponent;
    int    kappa,l; // l is redundant but convenient for r^l
    double itsNormalization;
};

//
//  Derived from the LargeBF P(r)=r^l*exp(-e*r^2) as:
//
//                                   { -2e*r^(l+1)*exp(-e*r^2), kappa<0
//     Q(r)=(d/dr+(1+kappa)/r)P(r) = {
//                                   { ((2l+1)/r-2er)*r^l*exp(-e*r^2), kappa>0
//
class Small_BasisFunction
    : public virtual ::TBasisFunction<double>
{
public:
    Small_BasisFunction();
    Small_BasisFunction(const Large_BasisFunction*,double norm);

    virtual bool   operator==(const ::BasisFunction&) const;

    virtual std::ostream&       Write(std::ostream&) const;
    virtual std::istream&       Read (std::istream&)      ;
    virtual ::BasisFunction* Clone(        ) const;

    virtual double operator()(const Vec3&) const;
    virtual Vec3   Gradient  (const Vec3&) const;
private:
    const Large_BasisFunction* Pr;
    double itsNormalization;
};


}} //namespace

#endif //_Atom_kappa_Gaussian_BF_H_
