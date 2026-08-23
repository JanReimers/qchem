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
// i.e. the fit coefficients ARE the field's values at the mesh points, and the division is exact rather
// than merely convenient.  That identity is why this fitter samples and stores, where its orthoNORMAL
// sibling samples and forward-FFTs: same shape, and the diagonal it divides by is the one that sibling
// asserts to be all ones (OrthoNormalFunctionFitter.C).  Sampling is done BY THE BASIS -- it owns the
// points -- so nothing here knows where they are or what kind of mesh they came from.
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

    //! The fit: the BASIS samples the field at its own points, and the coefficients are those values --
    //! \f$c_g=\langle\delta_g|f\rangle/w_g\f$, with the \f$w_g\f$ cancelling exactly (see the file header).
    //! Written as the sample rather than as multiply-then-divide precisely BECAUSE it cancels: forming
    //! \f$w_g f_g\f$ and dividing it back would move the last bit of every coefficient for nothing.
    virtual void DoFit(const ProjectedScalar_R& ps) override
    {
        itsC=itsFitBasis->Sample(*ps.GetScalarFunction());
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
