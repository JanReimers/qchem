// File: BSpline/IE_Primatives.H get all calculation of primative integrals in one place.
#ifndef _BSpline_IE_Primatives_H_
#define _BSpline_IE_Primatives_H_

#include "IE.H"

namespace BSpline
{
template <size_t K> class IE_Primatives
    : public virtual Primative_Overlap  <double,K>
    , public virtual Primative_Grad2    <double,K>
    , public virtual Primative_Nuclear  <double,K>
    , public virtual Primative_Repulsion<double,K>
    , public virtual Primative_Charge   <double,K>
{
    typedef bspline::Spline<double, K> spline_t;
protected:
    virtual double Overlap  (const spline_t& a , const spline_t& b,size_t l_total) const;
    virtual double Grad2    (const spline_t& a , const spline_t& b,size_t la, size_t lb) const;
    virtual double Nuclear  (const spline_t& a , const spline_t& b,size_t l_total) const;
    virtual double Repulsion(const spline_t& a , const spline_t& b,size_t la,size_t lc) const;
    virtual double Charge   (const spline_t& a ,                   size_t l) const;
};

}

#endif // _BSpline_IE_Primatives_H_
   