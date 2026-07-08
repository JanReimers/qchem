// File: Fitting/Internal/OrthoFunctionFitter.C  Orthonormal (plane-wave / G-space) density fitter.
//
// The orthonormal sibling of the molecular ConstrainedFF: on an orthonormal {G} fit basis the density's
// rho-tilde IS the fit (no metric solve), so this implements ONLY the minimal CORE FunctionFitter_Density
// face -- no self-energy / charge constraint / rescale / real-space value (the non-ortho refinement an AO
// fit needs).  DoFit RECEIVES the density's pre-computed rho-tilde (a ProjectedDensity_G) and Repulsion
// delegates the FFT-free G-space Poisson solve to the orbital basis (Band_FT_IBS).
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
import qchem.BasisSet.Band_FT_IBS;            // the orbital assembly bridge (MakePotential) for the XC matrix
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
    //! Repulsion3C tensor (D contracted to V_H, kernel baked) and assembled by PW_Hartree via MakePotential --
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
    : public virtual GriddedScalarFitter
{
public:
    typedef std::shared_ptr<const BasisSet::cFIT_SF_ABS> fbs_t;
    explicit OrthoScalarFitter(const fbs_t& fbs) : itsFitBasis(fbs) {}

    //! GriddedScalarFitter: the fitter's own FFT quadrature grid engine (reached from the neutral fit face we
    //! hold).  The XC term borrows this ONE grid rather than cross-casting the fit basis a second time (#7).
    virtual const BasisSet::G_FieldEvaluator& Grid() const override {return FitGrid();}

    //! The "fit": batch-sample v_xc(r) on the fit basis's own FFT grid, then forward-transform.  Orthonormal
    //! exactness => the projection IS the fit (no metric solve); the field's FFT fast path is its own business.
    virtual void DoFit(const ProjectedScalar_R& ps) override
    {
        const BasisSet::G_FieldEvaluator& ge=FitGrid();
        rvec_t vals=(*ps.GetScalarFunction())(ge.GridPoints());   // batch-sample v_xc on the fit grid
        itsVt =ge.ForwardFFT(vals);                               // full /Npts G grid (for the assembly)
        itsMap=ge.FieldCoeffs(itsVt);                             // fit-basis coefficients (for op(r) plotting)
    }

    //! XC matrix <i|v_xc|j> = V-tilde(m_i-m_j): the ORBITAL basis assembles over ITS {G}, looking each
    //! reciprocal-index difference up in OUR fit-grid coefficients -- so the fit grid may be denser than (or
    //! offset from) the orbital's.  The orbital's assembly is its Band_FT_IBS::MakePotential bridge (the
    //! orbital-specific step); the fit grid's GridCoeff lookup is the shared G_FieldEvaluator grid engine.
    virtual hmat_t<dcmplx> Overlap(const robs_t<dcmplx>* bs) const override
    {
        // The orbital basis is CALLER-supplied and carries no compile-time guarantee of being plane-wave, so
        // this is a genuine "is it?" cross-cast: a reference-cast THROWS std::bad_cast (not release-mode UB)
        // for any future non-PW complex orbital basis.  (Contrast the fitter's own itsFitBasis casts, which
        // its isOrtho() contract guarantees.)  Ties to the item-C dynamic_cast survey.
        const BasisSet::Band_FT_IBS&      orb=dynamic_cast<const BasisSet::Band_FT_IBS&>(*bs);   // the assembly bridge
        const BasisSet::G_FieldEvaluator& fit=FitGrid();                                         // the fit grid engine
        return orb.MakePotential([&](const ivec3_t& dm)->dcmplx {return fit.GridCoeff(itsVt, dm);});
    }

    virtual void ReScale(double factor) override
    {
        itsVt *= factor;
        for (auto& kv : itsMap) kv.second *= factor;
    }

    //! ScalarFunction (core): the fitted POTENTIAL v_xc,fit(r) = Re Σ_G V-tilde(G) e^{iG·r} over the fit {G},
    //! via the basis's G_FieldEvaluator -- what the GUI plots; the seed of the v_xc - v_xc,fit fit-residual.
    virtual double  operator()(const rvec3_t& r) const override {return FitGrid().EvalField(itsMap, r);}
    virtual rvec3_t Gradient  (const rvec3_t& r) const override {return FitGrid().EvalFieldGradient(itsMap, r);}

    virtual std::ostream& Write(std::ostream& os) const override
        {return os << "OrthoScalarFitter (plane-wave grid quadrature)" << std::endl;}

private:
    //! The fit basis's FFT grid engine (the DIP seam), reached from the neutral fit face we hold.
    const BasisSet::G_FieldEvaluator& FitGrid() const
    {
        auto* ge=dynamic_cast<const BasisSet::G_FieldEvaluator*>(itsFitBasis.get());
        assert(ge && "OrthoScalarFitter: the {G} fit basis must provide the G_FieldEvaluator grid engine");
        return *ge;
    }
    fbs_t      itsFitBasis;   //!< the {G} fit basis -- owns the (possibly denser) quadrature grid
    cvec_t     itsVt;         //!< v_xc forward-FFT'd on the fit grid (full raster grid), for the operator assembly
    ΔG_Map     itsMap;        //!< the fit-basis coefficients, for op(r)/Gradient (GUI / fit-residual)
};

} //namespace
