// File: Fitting/Internal/OrthoNormalFunctionFitter.C  Scalar fitter for an ORTHONORMAL fit basis.
//
// THE METRIC IT ASSUMES IS IN ITS NAME (user, 2026-08-22).  isOrtho() asks ORTHOGONAL -- metric DIAGONAL,
// so no linear SOLVE -- and that admits two cases: orthoNORMAL (metric = I), where the projection IS the
// fit outright, and merely orthogonal (metric = diag(<f_a|f_a>)), where the fit is that projection DIVIDED
// by the diagonal.  A plane-wave {G} basis over the cell is the first; a delta basis, with
// <delta_g|delta_g> = w_g, is the second.  Dropping the denominator is right for one and an error of a
// factor w_g for the other, which is why this fitter says which it is:
//
//   IT ASSUMES OverlapDiagonal() == {1,1,1,...} AND NEVER CALLS IT.
//
// That is not laziness -- it is the specialisation earning its file.  For a raster fit basis the diagonal
// is npts ones nobody should materialise, and the forward FFT already delivers the coefficients directly.
// A general orthogonal fitter (the delta case) is its sibling: same shape, one extra elementwise divide.
module;
#include <cassert>
#include <memory>
#include <ostream>
export module qchem.Fitting.Internal.OrthoNormalFunctionFitter;
export import qchem.Fitting.FunctionFitter;  // FunctionFitter_Scalar + FitContraction<dcmplx,dcmplx>, ProjectedScalar_R
import qchem.Fitting.Types;                   // robs_t<dcmplx>
import qchem.BasisSet.Orbital_DFT_IBS;                // cFIT_SF_ABS (the held fit basis)
import qchem.BasisSet.G_FieldEvaluator;       // the DIP seam: inverse-transform itsMap to real space (op(r))
import qchem.Blaze;                           // hmat_t<dcmplx>

export namespace qchem::Fitting
{

//! \brief Scalar (overlap-metric) fitter on an orthonormal (plane-wave, G-space) fit basis -- the minimal
//! CORE face.  The XC sibling of OrthoFunctionFitter: DoFit SAMPLES the real field v_xc(r) on the fit basis's
//! OWN FFT grid (the G_FieldEvaluator grid engine, at the fit basis's -- possibly denser -- Ecut) and forward-
//! transforms it; Overlap has the ORBITAL basis assemble <i|v_xc|j> = V-tilde(m_i-m_j), looking each difference
//! up in the fit-grid coefficients.  This is what makes the fit grid come from the FIT basis (not the orbital):
//! relCutoff / the functional's GridCutoffFactor now actually control the XC quadrature.  Created through
//! Factory(cFIT_SF_ABS).
class OrthoNormalScalarFitter
    : public virtual FunctionFitter_Scalar
    , public virtual FitContraction<dcmplx,dcmplx>   // a raster fit contracts COMPLEX blocks on a COMPLEX fit axis
    , public virtual ScalarFunction<double>   // ...and its fit IS an evaluatable field (inverse transform)
{
public:
    typedef std::shared_ptr<const BasisSet::cFIT_SF_ABS> fbs_t;
    explicit OrthoNormalScalarFitter(const fbs_t& fbs) : itsFitBasis(fbs) {}

    //! The "fit": batch-sample v_xc(r) on the fit basis's own RASTER, then forward-transform.  Orthonormal
    //! exactness => the projection IS the fit (no metric solve).
    //!
    //! \note It does NOT go through \c FIT_SF_ABS::Overlap, and the reason is not bit-preservation
    //! (2026-08-23).  That face projects onto the fit basis's \f$\{G\}\f$ BALL -- one coefficient per fit
    //! FUNCTION, zero outside -- whereas \c Overlap(bs) below looks \f$\tilde v\f$ up at ORBITAL index
    //! differences \f$m_i-m_j\f$, which run to twice the orbital ball and are answered by the RASTER's
    //! wrap-around (\c GridCoeff).  Fitting through the ball would silently delete most of \f$v_{xc}\f$.
    //! So this fitter's "fit basis" for assembly purposes is the raster's \f$\{G\}\f$, and it says so by
    //! asking the raster face for both halves.  (\c G_FieldEvaluator::ProjectField carries the same
    //! truncation-vs-aliasing warning in the abstract; this is where it bites.)
    virtual void DoFit(const ProjectedScalar_R& ps) override
    {
        const BasisSet::G_RasterTransform& ge=FitRaster();
        itsVt =ge.ForwardFFT(ge.Sample(*ps.GetScalarFunction())); // full /Npts G grid (for the assembly)
        itsMap=ge.FieldCoeffs(itsVt);                             // fit-basis coefficients (for op(r) plotting)
    }

    //! XC matrix <i|v_xc|j> = V-tilde(m_i-m_j): the ORBITAL basis assembles over ITS {G}, looking each
    //! reciprocal-index difference up in OUR fit-grid coefficients -- so the fit grid may be denser than (or
    //! offset from) the orbital's.  The orbital's assembly is its applyAdjoint realization (the
    //! orbital-specific step); the fit grid's GridCoeff lookup is the shared G_FieldEvaluator grid engine.
    virtual hmat_t<dcmplx> Overlap(const BasisSet::Orbital_DFT_IBS<dcmplx,dcmplx>& orb) const override
    {
        // ASK THE FIT BASIS for <i|f_a|j> (2026-08-23): FIT_SF_ABS::Overlap3C's default delegates straight
        // back to orb.Overlap3C(*this), so this is the same cached tensor this line used to fetch by hand.
        const BasisSet::G_RasterTransform& fit=FitRaster();                           // the coefficient lookup
        // <i|v_xc|j> = Σ_k v_xc-tilde(G_k) <i|e^{iG_k}|j> -- the BACKWARD contraction of the OVERLAP 3-centre
        // tensor over OUR fit basis (which carries the fit {G}/grid), so the KS matrix is integrated back on the
        // SAME grid the density was collocated on (doc/GPWPlan §0e step 2).  No grid-less MakeOverlap(field).
        // The CONTRACTION FACE takes the DFT block since 2026-08-24, so "is it DFT-capable?" is answered by
        // the compiler at the call site and there is no cast in this fitter at all.
        return ContractAdjoint(itsFitBasis->Overlap3C(orb),
                                     [&](const ivec3_t& dm)->dcmplx {return fit.GridCoeff(itsVt, dm);});
    }

    virtual void ReScale(double factor) override
    {
        itsVt *= factor;
        for (auto& kv : itsMap) kv.second *= factor;
    }

    //! ScalarFunction (core): the fitted POTENTIAL v_xc,fit(r) = Re Σ_G V-tilde(G) e^{iG·r} over the fit {G},
    //! via the basis's G_FieldEvaluator -- what the GUI plots; the seed of the v_xc - v_xc,fit fit-residual.
    virtual double  operator()(const rvec3_t& r) const override {return FieldEval().EvalField(itsMap, r);}
    virtual rvec3_t Gradient  (const rvec3_t& r) const override {return FieldEval().EvalFieldGradient(itsMap, r);}

    virtual std::ostream& Write(std::ostream& os) const override
        {return os << "OrthoNormalScalarFitter (orthonormal {G}: the projection IS the fit)" << std::endl;}

private:
    //! The fit basis's RASTER TRANSFORMS -- the FFT half, which only a raster-backed basis can answer.
    const BasisSet::G_RasterTransform& FitRaster() const
    {
        auto* ge=dynamic_cast<const BasisSet::G_RasterTransform*>(itsFitBasis.get());
        assert(ge && "OrthoNormalScalarFitter: the {G} fit basis must provide the G_RasterTransform FFT pair");
        return *ge;
    }
    //! The evaluate face for op(r)/Gradient -- the same engine, the narrower ask.
    const BasisSet::G_FieldEvaluator& FieldEval() const
    {
        auto* fe=dynamic_cast<const BasisSet::G_FieldEvaluator*>(itsFitBasis.get());
        assert(fe && "OrthoNormalScalarFitter: the {G} fit basis must provide G_FieldEvaluator to evaluate v_xc,fit(r)");
        return *fe;
    }
    fbs_t      itsFitBasis;   //!< the {G} fit basis -- owns the (possibly denser) quadrature grid
    cvec_t     itsVt;         //!< v_xc forward-FFT'd on the fit grid (full raster grid), for the operator assembly
    ΔG_Map     itsMap;        //!< the fit-basis coefficients, for op(r)/Gradient (GUI / fit-residual)
};

} //namespace
