// File: Atom/l/BSpline_BF.H  B-Spline basis function.
#ifndef _Atom_l_BSpline_BF_H_
#define _Atom_l_BSpline_BF_H_

#include <BasisFunction.H>
#include <bspline/Core.h>

namespace Atoml
{
namespace BSpline
{
    
template <size_t K> class BasisFunction
    : public virtual ::TBasisFunction<double>
{
protected:
    typedef bspline::Spline<double, K> spline_t;
public:
    // BasisFunction();
    BasisFunction(const spline_t&, int l, double norm);
    
    virtual std::ostream&    Write(std::ostream&) const;
    virtual ::BasisFunction* Clone(        ) const;

    virtual double operator()(const Vec3&) const;
    virtual Vec3   Gradient  (const Vec3&) const;

private:
    spline_t itsSpline;
    bspline::Spline<double, K-1> itsDxSpline;
    int      itsL;
    double   itsNormalization;
};


}} //namespace Atom_l::Slater

#endif //_Atom_l_BSpline_BF_H_
