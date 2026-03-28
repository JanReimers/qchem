// File: BasisSet/Atom/BSpline/NR/IBS_Evaluator.C
module;
#include <bspline/Core.h>
#include <valarray>
#include <vector>
#include <memory>
#include <iosfwd>
export module BasisSet.Atom.BSpline.NR.IBS_Evaluator;
import qchem.BasisSet.Atom.IBS_Evaluator;
import qchem.Basisset.Atom.BSpline.GLQuadrature;

class BSplineTests;
export template <size_t K> class BSpline_IBS : public IBS_Evaluator
{
    typedef bspline::Spline<double, K> spline_t;
public: 
 
    BSpline_IBS(size_t Ngrid, double rmin, double rmax, int _l, const is_t& _mls);
    BSpline_IBS(size_t Ngrid, double rmin, double rmax, int _l) : BSpline_IBS(Ngrid,rmin,rmax,_l,{}) {};
    virtual void Register(Grouper*); //Set up unique spline or exponent indexes.
    
    virtual std::ostream& Write   (std::ostream&) const;
    virtual size_t maxSpan() const {return l<=K ? K-l : 0;}  //assume no overlap for indices separated by > maxSpan

    virtual rsmat_t Overlap  () const;
    virtual rsmat_t Grad2    () const;
    virtual rsmat_t Inv_r1   () const;
    virtual rsmat_t Inv_r2   () const;
    virtual rsmat_t Repulsion() const;
    virtual  rvec_t Charge   () const;
    virtual ds_t   Norm     () const {return ns;}
    virtual rmat_t XRepulsion(const Fit_IBS&) const;
    virtual rmat_t XKinetic  (const Orbital_RKBS_IBS<double>*) const;

    virtual dERI3  Overlap  (const Fit_IBS&) const; //3 center
    virtual dERI3  Repulsion(const Fit_IBS&) const; //3 center

    virtual rvec_t     operator() (const rvec3_t&) const;
    virtual rvec3vec_t Gradient   (const rvec3_t&) const;

    const spline_t& operator[](int index) const {return splines[index];}

protected:
    friend class BSplineTests;
    std::vector<double> MakeLogKnots(size_t NGrid, double rmin, double rmax);
    ds_t norms() const; //assumes es,l are already initialized

    double rmin,rmax; //This might be needed for creating fit basis sets.
    std::vector<spline_t> splines;
    std::unique_ptr<GLCache> itsGL;
};