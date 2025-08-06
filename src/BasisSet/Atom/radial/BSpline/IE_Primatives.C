// File: BSpline/IE_Primatives.C get all calculation of primative integrals in one place.
module;
#include <bspline/Core.h>
export module qchem.Basisset.Atom.radial.BSpline.IE_Primatives;
import qchem.Basisset.Atom.radial.BSpline.IE;

export namespace BSpline
{
template <size_t K> class IE_Primatives_Imp
    : public virtual IE_Primatives<K>
    , public virtual Primative_Overlap  <double,K>
    , public virtual Primative_Grad2    <double,K>
    , public virtual Primative_Inv_r1   <double,K>
    , public virtual Primative_Inv_r2   <double,K>
    , public virtual Primative_Repulsion<double,K>
    , public virtual Primative_Charge   <double,K>
{
    typedef bspline::Spline<double, K> spline_t;
protected:
    virtual double Overlap  (const spline_t& a , const spline_t& b,size_t l_total) const;
    virtual double Grad2    (const spline_t& a , const spline_t& b,size_t la, size_t lb) const;
    virtual double Inv_r1   (const spline_t& a , const spline_t& b,size_t l_total) const; //! <a|1/r^1|b>
    virtual double Inv_r2   (const spline_t& a , const spline_t& b,size_t l_total) const; //! <a|1/r^2|b>
    virtual double Repulsion(const spline_t& a , const spline_t& b,size_t la,size_t lc) const;
    virtual double Charge   (const spline_t& a ,                   size_t l) const;
};

}
   