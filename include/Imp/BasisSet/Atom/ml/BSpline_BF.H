// File: Atom/ml/BSpline_BF.H  B-Spline basis function.
#ifndef _Atom_ml_BSpline_BF_H_
#define _Atom_ml_BSpline_BF_H_

#include <bspline/Core.h>
#include <iosfwd>
import qchem.BasisFunction;

namespace Atom_ml
{
namespace BSpline
{
    
template <size_t K> class BasisFunction
    : public Atoml::BSpline::BasisFunction<K>
{
    typedef Atoml::BSpline::BasisFunction<K> Base;
    typedef typename Base::spline_t spline_t;
public:
    // BasisFunction();
    BasisFunction(const spline_t&, int l, int ml, double norm);
    
    virtual ::Real_BF* Clone(        ) const;

private:
    int      ml;
};


}} //namespace Atom_ml::Slater

#endif //_Atom_ml_BSpline_BF_H_
