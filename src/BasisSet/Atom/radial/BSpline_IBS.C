// File: BasisSet/Atom/radial/BSpline_IBS.C
module;
#include <bspline/Core.h>
#include <valarray>
#include <vector>
#include <memory>
export module BasisSet.Atom.BSpline_IBS;
import BasisSet.Atom.IBS_Evaluator;
import qchem.Basisset.Atom.radial.BSpline.GLQuadrature;


export template <size_t K> class BSpline_IBS : public virtual IBS_Evaluator
{
    typedef bspline::Spline<double, K> spline_t;
public: 
 
    BSpline_IBS(size_t Ngrid, double rmin, double rmax, int _l, const is_t& _mls);
    virtual size_t size() const {return splines.size();}

    virtual omls_t Overlap  () const;
    virtual omls_t Grad2    () const;
    virtual omls_t Inv_r1   () const;
    virtual omls_t Inv_r2   () const;
    virtual omls_t Repulsion() const;
    virtual omlv_t Charge   () const;

    // virtual ERI3   Overlap  (const BSpline_IBS&) const; //3 center
    // virtual ERI3   Repulsion(const BSpline_IBS&) const; //3 center

    virtual Vec     operator() (const RVec3&) const;
    virtual Vec3Vec Gradient   (const RVec3&) const;

private:
    std::vector<double> MakeLogKnots(size_t NGrid, double rmin, double rmax);
    ds_t norms() const; //assumes es,l are already initialized

    double rmin,rmax; //This might be needed for creating fit basis sets.
    std::vector<spline_t> splines;
    std::unique_ptr<GLCache> itsGL;
    int  l;
    is_t mls;
    ds_t ns;
};