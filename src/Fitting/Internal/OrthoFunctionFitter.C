// File: Fitting/Internal/OrthoFunctionFitter.C  Orthonormal (plane-wave / G-space) DENSITY fitter.
//
// The SCALAR (potential) sibling moved to OrthoNormalFunctionFitter.C 2026-08-22, where its name states
// the metric it assumes: it is orthoNORMAL-only, and the distinction became load-bearing when a merely
// ORTHOGONAL fit basis (delta, metric diag(w_g)) started answering isOrtho()==true.
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
import qchem.BasisSet.Orbital_DFT_IBS;                // cFIT_CD_ABS / cFIT_SF_ABS (the held fit bases)
import qchem.BasisSet.G_FieldEvaluator;       // the DIP seam: inverse-transform itsMap to real space (op(r))
import qchem.Blaze;                           // hmat_t<dcmplx>

export namespace qchem::Fitting
{

//! \brief Density fitter on an orthonormal (plane-wave, G-space) fit basis -- the minimal CORE face only.
class OrthoFunctionFitter
    : public virtual FunctionFitter_Density<dcmplx>
    , public virtual ScalarFunction<double>   // ...and ITS fit IS a field: the inverse transform over {G}
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

} //namespace
