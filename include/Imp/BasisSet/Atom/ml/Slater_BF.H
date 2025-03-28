// File: Atom/ml/Slater_BF.H  r^l exp(-ar)*Y_lm type basis function 
#ifndef _Atom_ml_Slater_BF_H_
#define _Atom_ml_Slater_BF_H_

#include <BasisFunction.H>

namespace Atom_ml
{
namespace Slater
{
    
class BasisFunction
    : public virtual ::TBasisFunction<double>
{
public:
    BasisFunction();
    BasisFunction(double ex, int n, int l, int m, double norm);
    
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


}} //namespace _Atom_ml_Slater_BF_H_

#endif //_SlaterBF_H_
