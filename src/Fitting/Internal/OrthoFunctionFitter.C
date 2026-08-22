// File: Fitting/Internal/OrthoFunctionFitter.C  Orthonormal (plane-wave / G-space) density fitter.
//
// The orthonormal sibling of the molecular ConstrainedFF: on an orthonormal {G} fit basis the density's
// rho-tilde IS the fit (no metric solve), so this implements ONLY the minimal CORE FunctionFitter_Density
// face -- no self-energy / charge constraint / rescale / real-space value (the non-ortho refinement an AO
// fit needs).  DoFit RECEIVES the density's pre-computed rho-tilde (a ProjectedDensity_G) and Repulsion
// delegates the FFT-free G-space Poisson solve to the orbital basis (Orbital_DFT_IBS<dcmplx>).
//
// It is created through the factory Factory(cFIT_CD_ABS) exactly as the AO fitter is created
// through Factory(rFIT_CD_ABS): the plane-wave Hartree term obtains it via the basis's
// CreateCDFitBasisSet, never assuming orbital==fit.  It HOLDS the fit basis (the tunable {G} grid) for the
// future denser-grid resampling; today that resampling is the identity, so the held basis is inert.
module;
#include <cassert>
#include <memory>
#include <ostream>
#include <stdexcept>
export module qchem.Fitting.Internal.OrthoFunctionFitter;
export import qchem.Fitting.FunctionFitter;  // FunctionFitter_Density/_Scalar<dcmplx>, ProjectedDensity/Scalar_G, ΔG_Map
import qchem.Fitting.Types;                   // robs_t<dcmplx>
import qchem.BasisSet.Fit_IBS;                // cFIT_CD_ABS / cFIT_SF_ABS (the held fit bases)
import qchem.BasisSet.G_FieldEvaluator;       // the DIP seam: inverse-transform itsMap to real space (op(r))
import qchem.BasisSet.Quadrature;             // BasisSet::Quadrature -- where the fit basis says it samples
import qchem.BasisSet.Orbital_DFT_IBS;            // the orbital assembly bridge (MakeOverlap) for the XC matrix
import qchem.Blaze;                           // hmat_t<dcmplx>

export namespace qchem::Fitting
{

//! \brief Density fitter on an orthonormal (plane-wave, G-space) fit basis -- the minimal CORE face only.
class OrthoFunctionFitter
    : public virtual FunctionFitter_Density<dcmplx>
{
public:
    typedef std::shared_ptr<const BasisSet::cFIT_CD_ABS> fbs_t;
    explicit OrthoFunctionFitter(const fbs_t& fbs) : itsFitBasis(fbs) {}

    //! The "fit": receive the density's G-space coefficients (rho-tilde).  Orthonormal exactness => nothing
    //! to solve, just store; the neutral argument is a ProjectedDensity_G (a sanctioned cross-cast).
    virtual void DoFit(const ProjectedDensity<dcmplx>& pd) override
    {
        auto g=dynamic_cast<const ProjectedDensity_G*>(&pd);
        assert(g && "OrthoFunctionFitter::DoFit requires a ProjectedDensity_G (G-space) projection");
        itsMap=g->Map();
    }

    //! Coulomb (Hartree) matrix: OBSOLETE.  The plane-wave Hartree matrix is now built from the density's
    //! Repulsion3C tensor (D contracted to V_H, kernel baked) and assembled by PW_Hartree via MakeOverlap --
    //! the density fitter is no longer on that path (it survives only for its evaluatable rho_fit(r) field).
    //! Nothing calls this; kept to satisfy the FunctionFitter_Density interface, throwing if ever reached.
    virtual hmat_t<dcmplx> Repulsion(const robs_t<dcmplx>*) const override
    {
        throw std::logic_error("OrthoFunctionFitter::Repulsion is obsolete: the plane-wave Hartree matrix comes "
                               "from the density's Repulsion3C tensor (PW_Hartree), not the density fitter.");
    }

    //! ScalarFunction (core): the fitted DENSITY rho_fit(r) = Re Σ_dm rho-tilde(dm) e^{i(B·dm)·r}.  We hold
    //! only the structure-neutral {G} fit face, so we DELEGATE the inverse transform to the basis's
    //! G_FieldEvaluator (DIP) -- a sanctioned "I want more" cross-cast to a capability, never a cast to
    //! concrete PW.  (What the GUI plots; the seed of the fit-residual rho - rho_fit diagnostic.)
    virtual double  operator()(const rvec3_t& r) const override {return FieldEval().EvalField(itsMap, r);}
    virtual rvec3_t Gradient  (const rvec3_t& r) const override {return FieldEval().EvalFieldGradient(itsMap, r);}

    virtual std::ostream& Write(std::ostream& os) const override
        {return os << "OrthoFunctionFitter (orthonormal G-space fit)" << std::endl;}

private:
    //! Reach the fit basis's inverse-transform capability (the DIP seam) from the neutral face we hold.
    const BasisSet::G_FieldEvaluator& FieldEval() const
    {
        auto* fe=dynamic_cast<const BasisSet::G_FieldEvaluator*>(itsFitBasis.get());
        assert(fe && "OrthoFunctionFitter: the {G} fit basis must provide G_FieldEvaluator to evaluate rho_fit(r)");
        return *fe;
    }
    fbs_t      itsFitBasis;   //!< the tunable {G} fit basis (the factory seam; inert until denser-grid resampling)
    ΔG_Map     itsMap;        //!< the fit = the density's rho-tilde (received in DoFit)
};

//! \brief Scalar (overlap-metric) fitter on an orthonormal (plane-wave, G-space) fit basis -- the minimal
//! CORE face.  The XC sibling of OrthoFunctionFitter: DoFit SAMPLES the real field v_xc(r) on the fit basis's
//! OWN FFT grid (the G_FieldEvaluator grid engine, at the fit basis's -- possibly denser -- Ecut) and forward-
//! transforms it; Overlap has the ORBITAL basis assemble <i|v_xc|j> = V-tilde(m_i-m_j), looking each difference
//! up in the fit-grid coefficients.  This is what makes the fit grid come from the FIT basis (not the orbital):
//! relCutoff / the functional's GridCutoffFactor now actually control the XC quadrature.  Created through
//! Factory(cFIT_SF_ABS).
class OrthoScalarFitter
    : public virtual FunctionFitter_Scalar<dcmplx>
{
public:
    typedef std::shared_ptr<const BasisSet::cFIT_SF_ABS> fbs_t;
    explicit OrthoScalarFitter(const fbs_t& fbs) : itsFitBasis(fbs) {}

    //! The "fit": batch-sample v_xc(r) at the fit basis's own quadrature POINTS, then project.  Orthonormal
    //! exactness => the projection IS the fit (no metric solve).  The two lines below are the two halves:
    //! WHERE to sample is the mesh's business, and the forward FFT is the raster's batch form of the
    //! per-\f$G\f$ weight-vector contraction (the field's own FFT fast path is its own business).
    virtual void DoFit(const ProjectedScalar_R& ps) override
    {
        const BasisSet::G_RasterTransform& ge=FitRaster();
        rvec_t vals=(*ps.GetScalarFunction())(FitMesh().Points());   // batch-sample v_xc on the fit grid
        itsVt =ge.ForwardFFT(vals);                               // full /Npts G grid (for the assembly)
        itsMap=ge.FieldCoeffs(itsVt);                             // fit-basis coefficients (for op(r) plotting)
    }

    //! XC matrix <i|v_xc|j> = V-tilde(m_i-m_j): the ORBITAL basis assembles over ITS {G}, looking each
    //! reciprocal-index difference up in OUR fit-grid coefficients -- so the fit grid may be denser than (or
    //! offset from) the orbital's.  The orbital's assembly is its applyAdjoint realization (the
    //! orbital-specific step); the fit grid's GridCoeff lookup is the shared G_FieldEvaluator grid engine.
    virtual hmat_t<dcmplx> Overlap(const robs_t<dcmplx>* bs) const override
    {
        // The orbital basis is CALLER-supplied and carries no compile-time guarantee of being plane-wave, so
        // this is a genuine "is it?" cross-cast: a reference-cast THROWS std::bad_cast (not release-mode UB)
        // for any future non-PW complex orbital basis.  (Contrast the fitter's own itsFitBasis casts, which
        // its isOrtho() contract guarantees.)  Ties to the item-C dynamic_cast survey.
        const BasisSet::Orbital_DFT_IBS<dcmplx>&      orb=dynamic_cast<const BasisSet::Orbital_DFT_IBS<dcmplx>&>(*bs);   // the assembly bridge
        const BasisSet::G_RasterTransform&            fit=FitRaster();                           // the coefficient lookup
        // <i|v_xc|j> = Σ_k v_xc-tilde(G_k) <i|e^{iG_k}|j> -- the BACKWARD contraction of the OVERLAP 3-centre
        // tensor over OUR fit basis (which carries the fit {G}/grid), so the KS matrix is integrated back on the
        // SAME grid the density was collocated on (doc/GPWPlan §0e step 2).  No grid-less MakeOverlap(field).
        return ContractAdjoint(orb.Overlap3C(*itsFitBasis),
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
        {return os << "OrthoScalarFitter (plane-wave grid quadrature)" << std::endl;}

private:
    //! The fit basis's QUADRATURE -- points+weights, the universal half (the DIP seam), reached from the
    //! neutral fit face we hold.  Guaranteed by Factory(cFIT_SF_ABS)'s construction-time contract, so
    //! asserts (not throws) suffice.
    const qcMesh::Mesh& FitMesh() const
    {
        auto* q=dynamic_cast<const BasisSet::Quadrature*>(itsFitBasis.get());
        assert(q && "OrthoScalarFitter: the {G} fit basis must carry its quadrature (BasisSet::Quadrature)");
        return q->Mesh();
    }
    //! The fit basis's RASTER TRANSFORMS -- the FFT half, which only a raster-backed basis can answer.
    const BasisSet::G_RasterTransform& FitRaster() const
    {
        auto* ge=dynamic_cast<const BasisSet::G_RasterTransform*>(itsFitBasis.get());
        assert(ge && "OrthoScalarFitter: the {G} fit basis must provide the G_RasterTransform FFT pair");
        return *ge;
    }
    //! The evaluate face for op(r)/Gradient -- the same engine, the narrower ask.
    const BasisSet::G_FieldEvaluator& FieldEval() const
    {
        auto* fe=dynamic_cast<const BasisSet::G_FieldEvaluator*>(itsFitBasis.get());
        assert(fe && "OrthoScalarFitter: the {G} fit basis must provide G_FieldEvaluator to evaluate v_xc,fit(r)");
        return *fe;
    }
    fbs_t      itsFitBasis;   //!< the {G} fit basis -- owns the (possibly denser) quadrature grid
    cvec_t     itsVt;         //!< v_xc forward-FFT'd on the fit grid (full raster grid), for the operator assembly
    ΔG_Map     itsMap;        //!< the fit-basis coefficients, for op(r)/Gradient (GUI / fit-residual)
};

} //namespace
