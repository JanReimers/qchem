// File: Atom/l/BSpline_BF.C  B-Spline basis function.

#include "Imp/BasisSet/Atom/l/BSpline_BF.H"
#include "Imp/Misc/IntPower.H"
#include <bspline/operators/Derivative.h>
#include "oml/vector3d.h"
#include <iostream>
#include <cassert>

namespace Atoml
{
namespace BSpline
{

//     template <size_t K> BasisFunction<K>::BasisFunction()
//     : itsSpline       ()
//     , itsL            (0)
//     , itsNormalization(0)
// {};

template <size_t K> BasisFunction<K>::BasisFunction(const spline_t& sp, int l, double norm)
    : itsSpline       (sp)
    , itsDxSpline     (transformSpline(bspline::operators::Dx<1>{},sp))
    , itsL            (l)
    , itsNormalization(norm)
{
    // bspline::operators::Derivative<1> ddx;
    // auto ddx1=bspline::operators::Dx<1>{};
    // auto dsdr=ddx1*(*itsSpline);
};

template <size_t K> bool BasisFunction<K>::operator==(const ::BasisFunction& bf) const
{
    const BasisFunction<K>& sbf = dynamic_cast<const BasisFunction<K>&>(bf);
    assert(&sbf);
    return (itsSpline==sbf.itsSpline) && (itsL==sbf.itsL);
}

template <size_t K> std::ostream& BasisFunction<K>::Write(std::ostream& os) const
{
    return os << "[" << itsSpline.front() << "," << itsSpline.back() << "] ";
}

template <size_t K> double BasisFunction<K>::operator()(const Vec3& r) const
{
    double mr=norm(r);
    if (mr<itsSpline.front() || mr>itsSpline.back()) return 0.0;
    return itsNormalization*itsSpline(mr);
}

template <size_t K> typename BasisFunction<K>::Vec3 BasisFunction<K>::Gradient(const Vec3& r) const
{
    Vec3 ret(0,0,0);
    if (r==ret) return ret; //Cusp at the origin so grad is undefined.
    double mr=norm(r);
    assert(mr>0);
    if (mr<itsSpline.front() || mr>itsSpline.back()) return ret;
    Vec3 rhat=r/mr;
    return rhat*itsNormalization*itsDxSpline(mr);
}

template <size_t K> ::BasisFunction* BasisFunction<K>::Clone() const
{
    return new  BasisFunction<K>(*this);
}

template class BasisFunction<6>;

}} //namespace
