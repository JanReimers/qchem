// File: Fitting/Internal/DeltaFunctionFitter.C  Scalar fitter on a DELTA (identity) fit basis.
//
// The third scalar fitter, and the one that makes the delta representation Liskov-substitutable with the
// other two: the XC term now obtains a fitter from the same Factory, calls the same DoFit/Overlap, and
// never learns which representation it got.
//
// WHAT IT IS.  A delta basis is ORTHOGONAL but not orthoNORMAL -- <delta_g|delta_h> = w_g delta_gh -- so
// the fit is the projection divided by the metric diagonal:
//
//     c_g = <delta_g|f> / <delta_g|delta_g> = w_g f(r_g) / w_g = f(r_g),
//
// i.e. the fit coefficients ARE the field's values at the mesh points -- but that is a FACT about this
// representation, not the shape of the code: the code is "projection, then my metric", the same two steps
// the Gaussian fitter takes with S^-1 and the orthoNORMAL sibling takes with nothing (the diagonal it
// asserts to be all ones -- OrthoNormalFunctionFitter.C).  Both halves are the BASIS's: it owns its
// functions and its metric, so nothing here knows where its points are or what kind of mesh they came from.
//
// AND WHY THE CONTRACTION IS THE BASIS'S.  Overlap is sum_g c_g <chi_i|delta_g|chi_j>, which is the
// weighted quadrature Phi^dag diag(w c) Phi -- the basis owns the points, the weights AND the Phi table,
// so it performs it and this class only decides WHICH c to hand over.  Both block scalars are served
// (FitContraction<double> and <dcmplx>), which a mixed real/complex run needs and which is exactly what
// the single-scalar fitter face could not express before the FitContraction split.
module;
#include <cassert>
#include <memory>
#include <ostream>
#include <stdexcept>
export module qchem.Fitting.Internal.DeltaFunctionFitter;
export import qchem.Fitting.FunctionFitter;  // FunctionFitter_Scalar + FitContraction<U>, ProjectedScalar_R
import qchem.Fitting.Types;                   // robs_t<U>
import qchem.BasisSet.Fit_IBS;                // cFIT_SF_Delta (the held fit basis)
import qchem.Blaze;                           // hmat_t<U>

export namespace qchem::Fitting
{

//! \brief Scalar (overlap-metric) fitter on a \f$\delta\f$ fit basis -- the projection IS the fit, up to
//! the basis's own metric diagonal.  Serves BOTH orbital-block scalars, as a mixed run requires.
class DeltaScalarFitter
    : public virtual FunctionFitter_Scalar
    , public virtual FitContraction<double>   // a real TRIM block's quadrature runs in real arithmetic
    , public virtual FitContraction<dcmplx>   // ...its complex siblings' in complex
{
public:
    typedef std::shared_ptr<const BasisSet::cFIT_SF_Delta> fbs_t;
    explicit DeltaScalarFitter(const fbs_t& fbs) : itsFitBasis(fbs)
    {
        assert(itsFitBasis && itsFitBasis->isOrtho() &&
               "DeltaScalarFitter: a delta basis is orthogonal (diagonal metric) by construction");
    }

    //! The fit, in the ONE shape every scalar fitter has: take the basis's projection and apply my metric.
    //! Here the metric is diagonal, so \f$c_a=\langle f_a|f\rangle/\langle f_a|f_a\rangle\f$ -- which for
    //! \f$\delta\f$ is \f$w_g f(r_g)/w_g\f$.  \c BasisSet::OrthogonalFit is that one line, shared with the
    //! only other consumer of it (a matrix-free density expressing itself over a \f$\delta\f$ basis).
    //! \note This is where \c OverlapDiagonal() finally gets a production consumer.  The division does NOT
    //! cancel to the last bit -- \f$\mathrm{fl}(w_gf_g)/w_g\ne f_g\f$ in general -- so the coefficients can
    //! move by an ulp against the old direct sample.  Measured on the pinned Si two-route gate: no printed
    //! digit and no iteration count moved (doc/CleanupCandidates.md R1.0 increment 3).
    virtual void DoFit(const ProjectedScalar_R& ps) override
    {
        itsC=BasisSet::OrthogonalFit(*itsFitBasis, *ps.GetScalarFunction());
    }

    //! \copydoc Fitting::FitContraction::Overlap
    //! \f$\sum_g c_g\langle\chi_i|\delta_g|\chi_j\rangle=\Phi^\dagger\mathrm{diag}(w\,c)\Phi\f$ -- the
    //! basis performs it (it owns the weights and the table); this only supplies \f$c\f$.
    virtual hmat_t<double> Overlap(const robs_t<double>* bs) const override
        {assert(bs); return itsFitBasis->Quadrature(*bs, itsC);}
    virtual hmat_t<dcmplx> Overlap(const robs_t<dcmplx>* bs) const override
        {assert(bs); return itsFitBasis->Quadrature(*bs, itsC);}

    virtual void ReScale(double factor) override {itsC *= factor;}

    //! \warning NO evaluatable fitted field, by design.  A \f$\delta\f$ expansion has no value BETWEEN its
    //! points -- \f$\sum_g c_g\delta_g(r)\f$ is a distribution -- so this fitter does not derive the
    //! \c ScalarFunction capability its two siblings do, and nothing asks it to.  (The AO fitter evaluates
    //! its Gaussians, the plane-wave fitter inverse-transforms its coefficients; both are genuine fields.)
    virtual std::ostream& Write(std::ostream& os) const override
        {return os << "DeltaScalarFitter (delta basis: the fit coefficients ARE the sampled values)" << std::endl;}

private:
    fbs_t  itsFitBasis;   //!< the δ basis -- the points, the weights, the Φ tables and the quadrature
    rvec_t itsC;          //!< the fit = the field's values at its points (see DoFit)
};

} //namespace
