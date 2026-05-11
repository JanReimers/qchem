// File: BasisSet1/Atom/Evaluators/BSpline/IBS_Evaluator.C
module;
#include <bspline/Core.h>
#include <vector>
#include <memory>
#include <iosfwd>
export module qchem.BasisSet1.Atom.Evaluators.BSpline.IBS;
import qchem.BasisSet1.Atom.Evaluators.IBS;
import qchem.BasisSet1.Atom.Evaluators.BSpline.Internal.GLQuadrature;
import qchem.Symmetry.Yl;
//
//  This version is for phi(r) = sum(Bi(r),i)
// 
export template <size_t K> class BSpline_IBS_Evaluator : public IBS_Evaluator
{
    typedef bspline::Spline<double, K> spline_t;
public: 
    BSpline_IBS_Evaluator(size_t Ngrid, double rmin, double rmax, const Irrep_QNs::sym_t& ylm);
    BSpline_IBS_Evaluator(const BSpline_IBS_Evaluator& b) : BSpline_IBS_Evaluator(b.splines.size(),b.rmin,b.rmax,Irrep_QNs::sym_t(new Yl_Sym(0))) {};
    BSpline_IBS_Evaluator Rescale(double scale_factor) const
    {
        return BSpline_IBS_Evaluator(*this);
    }    
    virtual void Register(Grouper*); //Set up unique spline or exponent indexes.
    
    virtual std::ostream& Write   (std::ostream&) const;
    virtual size_t maxSpan() const {return l<=K ? K-l : 0;}  //assume no overlap for indices separated by > maxSpan

    virtual rsmat_t Overlap  () const;
    virtual rsmat_t Grad2    () const;
    virtual rsmat_t Inv_r1   () const;
    virtual rsmat_t Inv_r2   () const;
    virtual rsmat_t Repulsion() const;
    virtual  rvec_t Charge   () const;
    virtual  rvec_t Norm     () const {return ns;}
    virtual rmat_t XRepulsion(const IBS_Evaluator&) const;
    virtual rmat_t XKinetic  (const IBS_Evaluator*) const;

    virtual dERI3  Overlap  (const IBS_Evaluator&) const; //3 center
    virtual dERI3  Repulsion(const IBS_Evaluator&) const; //3 center

    virtual rvec_t     operator() (const rvec3_t&) const;
    virtual rvec3vec_t Gradient   (const rvec3_t&) const;

    virtual std::string RadialID() const;
    virtual std::string Name    () const;

    const spline_t& operator[](int index) const {return splines[index];}

protected:
    std::vector<double> MakeLogKnots(size_t NGrid, double rmin, double rmax);
    rvec_t norms() const; //assumes es,l are already initialized

    double rmin,rmax; //This might be needed for creating fit basis sets.
    std::vector<double> knots;
    std::vector<spline_t> splines;
    std::unique_ptr<GLCache> itsGL;
};

export template <size_t K> class BSpline_r_IBS_Evaluator : public IBS_Evaluator
{
    typedef bspline::Spline<double, K> spline_t;
public: 
    BSpline_r_IBS_Evaluator(size_t Ngrid, double rmin, double rmax, const Irrep_QNs::sym_t& ylm);
    virtual void Register(Grouper*); //Set up unique spline or exponent indexes.
    
    virtual std::ostream& Write   (std::ostream&) const;
    virtual size_t maxSpan() const {return l<=K ? K-l : 0;}  //assume no overlap for indices separated by > maxSpan

    virtual rsmat_t Overlap  () const;
    virtual rsmat_t Grad2    () const;
    virtual rsmat_t Inv_r1   () const;
    virtual rsmat_t Inv_r2   () const;
    virtual rsmat_t Repulsion() const;
    virtual  rvec_t Charge   () const;
    virtual  rvec_t Norm     () const {return ns;}
    virtual rmat_t XRepulsion(const IBS_Evaluator&) const;
    virtual rmat_t XKinetic  (const IBS_Evaluator*) const;

    virtual dERI3  Overlap  (const IBS_Evaluator&) const; //3 center
    virtual dERI3  Repulsion(const IBS_Evaluator&) const; //3 center

    virtual rvec_t     operator() (const rvec3_t&) const;
    virtual rvec3vec_t Gradient   (const rvec3_t&) const;

    virtual std::string RadialID () const;
    virtual std::string Name    () const;

    const spline_t& operator[](int index) const {return splines[index];}

protected:
    std::vector<double> MakeLogKnots(size_t NGrid, double rmin, double rmax);
    rvec_t norms() const; //assumes es,l are already initialized

    double rmin,rmax; //This might be needed for creating fit basis sets.
    std::vector<double> knots;
    std::vector<spline_t> splines;
    std::unique_ptr<GLCache> itsGL;
};