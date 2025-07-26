// File: BSpline/IE_Primatives.C get all calculation of primative integrals in one place.
module;
#include <cassert>
#include <bspline/Core.h>
using namespace bspline::operators; 
using namespace bspline::integration; 
module qchem.Basisset.Atom.radial.BSpline.IE_Primatives;
import Common.Constants;
import qchem.Basisset.Atom.radial.BSpline.Rk;


namespace BSpline
{

template <size_t K> double IE_Primatives<K>::Repulsion(const spline_t& ab , const spline_t& c,size_t la,size_t lc) const
{    
    assert(false);
    return 0.0;
}

template <size_t K> double  IE_Primatives<K>::Overlap(const spline_t& a , const spline_t& b,size_t l_total) const
{
    return BilinearForm{X<2>{}}(a,b)*FourPi;
    // return BilinearForm{IdentityOperator{}}(a,b);
}

template <size_t K> double IE_Primatives<K>::Grad2(const spline_t& a , const spline_t& b,size_t la, size_t lb) const
{
    static const auto T = -X<2>{} * Dx<2>{} - 2 * X<1>{} * Dx<1>{};
    assert(la==lb);
    return BilinearForm{T}(a,b)*FourPi;
}

template <size_t K> double IE_Primatives<K>::Inv_r1(const spline_t& a , const spline_t& b,size_t l_total) const
{
    return BilinearForm{X<1>{}}(a,b)*FourPi; 
}
template <size_t K> double IE_Primatives<K>::Inv_r2(const spline_t& a , const spline_t& b,size_t l_total) const
{
    return BilinearForm{IdentityOperator{}}(a,b)*FourPi; 
}

template <size_t K> double IE_Primatives<K>::Charge(const spline_t& a , size_t l) const
{
    return LinearForm{X<2>{}}(a);
}

#define INSTANCEk(k) template class IE_Primatives<k>;
#include "../Instance.hpp"

} //namespace
