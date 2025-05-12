// File: BSpline/IE_Primatives.C get all calculation of primative integrals in one place.

#include "Imp/BasisSet/Atom/radial/BSpline/IE_Primatives.H"
#include "Imp/BasisSet/Atom/radial/BSpline/Integrals.H"
#include "Imp/BasisSet/Atom/radial/BSpline/Rk.H"

using namespace bspline::operators; 
using namespace bspline::integration; 

namespace BSpline
{

template <size_t K> double IE_Primatives<K>::Repulsion(const spline_t& ab , const spline_t& c,size_t la,size_t lc) const
{    
    assert(false);
    return 0.0;
}

template <size_t K> double  IE_Primatives<K>::Overlap(const spline_t& a , const spline_t& b,size_t l_total) const
{
    return BilinearForm{X<2>{}}(a,b)*4*Pi;
    // return BilinearForm{IdentityOperator{}}(a,b);
}

template <size_t K> double IE_Primatives<K>::Grad2(const spline_t& a , const spline_t& b,size_t la, size_t lb) const
{
    static const size_t L=0;
    static const auto T = -X<2>{} * Dx<2>{} - 2 * X<1>{} * Dx<1>{} + L*(L+1);
    assert(la==lb);
    assert(la==L);
    
    return BilinearForm{T}(a,b);
}

template <size_t K> double IE_Primatives<K>::Nuclear(const spline_t& a , const spline_t& b,size_t l_total) const
{
    return BilinearForm{X<1>{}}(a,b); //Already has 4*Pi
}

template <size_t K> double IE_Primatives<K>::Charge(const spline_t& a , size_t l) const
{
    return LinearForm{X<2>{}}(a);
}

template class IE_Primatives<6>;

} //namespace
