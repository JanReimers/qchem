// File: Atom/l/BSpline_BF.C  B-Spline basis function.
module;
#include <iostream>
#include <cassert>
#include <cmath>
#include <bspline/operators/Derivative.h>


module qchem.BasisSet.Atom.Internal.l.BSplineBS;

namespace Atoml
{
namespace BSpline
{

template <size_t K> BasisFunction<K>::BasisFunction(const spline_t& sp, int l, double norm)
    : itsSpline       (sp)
    , itsDxSpline     (transformSpline(bspline::operators::Dx<1>{},sp))
    , itsL            (l)
    , itsNormalization(norm)
{
};
template <size_t K> std::ostream& BasisFunction<K>::Write(std::ostream& os) const
{
    return os << "[" << itsSpline.front() << "," << itsSpline.back() << "] ";
}
template <size_t K> double BasisFunction<K>::operator()(const RVec3& r) const
{
    double mr=norm(r);
    if (mr<itsSpline.front() || mr>itsSpline.back()) return 0.0;
    return itsNormalization*itsSpline(mr);
}
template <size_t K> RVec3 BasisFunction<K>::Gradient(const RVec3& r) const
{
    RVec3 ret(0,0,0);
    if (r==ret) return ret; //Cusp at the origin so grad is undefined.
    double mr=norm(r);
    assert(mr>0);
    if (mr<itsSpline.front() || mr>itsSpline.back()) return ret;
    RVec3 rhat=r/mr;
    return rhat*itsNormalization*itsDxSpline(mr);
}

#define INSTANCEk(k) template class BasisFunction<k>;
#include "../../radial/BSpline/Instance.hpp"

}} //namespace
