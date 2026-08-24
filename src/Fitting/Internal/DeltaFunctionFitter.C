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
// (FitContraction<double,dcmplx> and <dcmplx,dcmplx>), which a mixed real/complex run needs and which is
// exactly what the single-scalar fitter face could not express before the FitContraction split.
module;
#include <cassert>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <type_traits>   // std::is_same_v (same-scalar face vs the mixed-run cross-cast)
export module qchem.Fitting.Internal.DeltaFunctionFitter;
export import qchem.Fitting.FunctionFitter;  // FunctionFitter_Scalar + FitContraction<U,TFit>, ProjectedScalar_R
import qchem.Fitting.Types;                   // robs_t<U>
import qchem.BasisSet.Orbital_DFT_IBS;                // cFIT_SF_ABS + Integrals_Overlap3C<U,TFit> (the held fit basis)
import qchem.Fitting.FitOperations;           // OrthogonalFit -- this library's projection/diagonal fit
import qchem.Blaze;                           // hmat_t<U>

export namespace qchem::Fitting
{

//! \brief Scalar (overlap-metric) fitter on a \f$\delta\f$ fit basis -- the projection IS the fit, up to
//! the basis's own metric diagonal.  Serves BOTH orbital-block scalars, as a mixed run requires.
class DeltaScalarFitter
    : public virtual FunctionFitter_Scalar
    , public virtual FitContraction<double,dcmplx>   // a real TRIM block's quadrature runs in real arithmetic
    , public virtual FitContraction<dcmplx,dcmplx>   // ...its complex siblings' in complex
{
public:
    typedef std::shared_ptr<const BasisSet::cFIT_SF_ABS> fbs_t;
    explicit DeltaScalarFitter(const fbs_t& fbs) : itsFitBasis(fbs)
    {
        assert(itsFitBasis && itsFitBasis->isOrtho() &&
               "DeltaScalarFitter: a delta basis is orthogonal (diagonal metric) by construction");
    }

    //! The fit, in the ONE shape every scalar fitter has: take the basis's projection and apply my metric.
    //! Here the metric is diagonal, so \f$c_a=\langle f_a|f\rangle/\langle f_a|f_a\rangle\f$ -- which for
    //! \f$\delta\f$ is \f$w_g f(r_g)/w_g\f$.  \c OrthogonalFit is that one line, shared with the
    //! only other consumer of it (a matrix-free density expressing itself over a \f$\delta\f$ basis).
    //! \note This is where \c OverlapDiagonal() finally gets a production consumer.  The division does NOT
    //! cancel to the last bit -- \f$\mathrm{fl}(w_gf_g)/w_g\ne f_g\f$ in general -- so the coefficients can
    //! move by an ulp against the old direct sample.  Measured on the pinned Si two-route gate: no printed
    //! digit and no iteration count moved (doc/CleanupCandidates.md R1.0 increment 3).
    virtual void DoFit(const ProjectedScalar_R& ps) override
    {
        itsC=OrthogonalFit(*itsFitBasis, *ps.GetScalarFunction());
    }

    //! \copydoc Fitting::FitContraction::Overlap
    //! \f$\langle\chi_i|\sum_g c_g\delta_g|\chi_j\rangle\f$: ask the basis for the 3-centre overlap over
    //! ITS functions and contract it with MY coefficients.  Structurally the same two calls the Gaussian
    //! fitter makes (\c Overlap3C then contract) -- which is what Liskov substitutability of the
    //! representations means in practice, rather than by analogy.
    virtual hmat_t<double> Overlap(const BasisSet::Orbital_DFT_IBS<double,dcmplx>& orb) const override {return Contract(orb);}
    virtual hmat_t<dcmplx> Overlap(const BasisSet::Orbital_DFT_IBS<dcmplx,dcmplx>& orb) const override {return Contract(orb);}

    virtual void ReScale(double factor) override {itsC *= factor;}

    //! \warning NO evaluatable fitted field, by design.  A \f$\delta\f$ expansion has no value BETWEEN its
    //! points -- \f$\sum_g c_g\delta_g(r)\f$ is a distribution -- so this fitter does not derive the
    //! \c ScalarFunction capability its two siblings do, and nothing asks it to.  (The AO fitter evaluates
    //! its Gaussians, the plane-wave fitter inverse-transforms its coefficients; both are genuine fields.)
    virtual std::ostream& Write(std::ostream& os) const override
        {return os << "DeltaScalarFitter (delta basis: the fit coefficients ARE the sampled values)" << std::endl;}

private:
    //! ONE body for both block scalars: the basis's \f$\langle\chi_i|\delta_g|\chi_j\rangle\f$, contracted
    //! against my coefficients.  A mixed run (3c-3) reaches it with either scalar.
    template <class U> hmat_t<U> Contract(const BasisSet::Orbital_DFT_IBS<U,dcmplx>& orb) const
    {
        // The 3-centre overlap is a DFT-TIER question, and since 2026-08-24 the CONTRACTION face says so
        // too: both overloads above take the DFT block, so the cross-cast that used to stand here belongs
        // to (and now lives with) the caller.  TFit is dcmplx throughout -- our fit basis is the periodic
        // one, and a real TRIM block against it is exactly the 3c-3 mixed case.
        // ...and MY fit basis's 3-centre face for THIS block's scalar.  For a Bloch block that face is a
        // (virtual) base of the fit face I hold and the cast never fails; for a real TRIM block it is the
        // genuine 3c-3 "can you serve a real block?" question, so by reference -- it throws rather than
        // being release-mode UB.  One line per site since 2026-08-24; it was a shared helper.
        const auto& face=dynamic_cast<const BasisSet::Integrals_Overlap3C<U,dcmplx>&>(*itsFitBasis);
        const Projector3<U>& o3=face.Overlap3C(orb);
        assert(o3.applyRawAdjoint && "DeltaScalarFitter: this fit basis must realise the 3-centre adjoint");
        return o3.applyRawAdjoint(itsC);
    }
    fbs_t  itsFitBasis;   //!< the δ basis -- my functions, their metric, and their 3-centre overlap
    rvec_t itsC;          //!< MY fit coefficients over that basis (see DoFit)
};

} //namespace
