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
export module qchem.Fitting.Internal.OrthoFunctionFitter;
export import qchem.Fitting.FunctionFitter;  // FunctionFitter_Density/_Scalar<dcmplx>, ProjectedDensity/Scalar_G, ΔG_Map
import qchem.Fitting.Types;                   // robs_t<dcmplx>
import qchem.BasisSet.Fit_IBS;                // cFIT_CD_ABS / cFIT_SF_ABS (the held fit bases)
import qchem.BasisSet.Band_FT_IBS;            // the reciprocal-space assembly the fits delegate to
import qchem.BasisSet.G_FieldEvaluator;       // the DIP seam: inverse-transform itsMap to real space (op(r))
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

    //! Coulomb (Hartree) matrix: delegate the FFT-free G-space Poisson solve to the orbital basis, which
    //! arrives as the common Orbital_1E_IBS base and is cast down to the reciprocal-space capability.
    virtual hmat_t<dcmplx> Repulsion(const robs_t<dcmplx>* bs) const override
    {
        auto pw=dynamic_cast<const BasisSet::Band_FT_IBS*>(bs);
        assert(pw && "OrthoFunctionFitter::Repulsion requires a Band_FT_IBS (plane-wave) basis");
        double Eh;   // discarded here: the Hartree ENERGY is the 1/2<rho|V_H> contraction the term does in
        return pw->Repulsion(itsMap, Eh);   // GetEnergy (DM_Contract), not this matrix-assembly by-product.
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
//! CORE face.  The XC sibling of OrthoFunctionFitter: DoFit RECEIVES the potential's
//! pre-computed V-tilde (a ProjectedScalar_G) and Overlap delegates the (kernel-free) assembly to the orbital
//! Band_FT_IBS.  Created through Factory(cFIT_SF_ABS); holds the {G} fit basis (the factory seam,
//! inert until the denser-{G} upgrade).
class OrthoScalarFitter
    : public virtual FunctionFitter_Scalar<dcmplx>
{
public:
    typedef std::shared_ptr<const BasisSet::cFIT_SF_ABS> fbs_t;
    explicit OrthoScalarFitter(const fbs_t& fbs) : itsFitBasis(fbs) {}

    //! The "fit": receive the potential's V-tilde (the term forward-FFT'd v_xc).  Orthonormal exactness =>
    //! nothing to solve, just store; the neutral argument is a ProjectedScalar_G (a sanctioned cross-cast).
    virtual void DoFit(const ProjectedScalar<dcmplx>& ps) override
    {
        auto g=dynamic_cast<const ProjectedScalar_G*>(&ps);
        assert(g && "OrthoScalarFitter::DoFit requires a ProjectedScalar_G (G-space) projection");
        itsMap=g->Map();
    }

    //! XC matrix <i|v_xc|j> = V-tilde(dm): assemble directly from the stored G-space coefficients (no kernel),
    //! delegating to the orbital basis's reciprocal-space assembly.
    virtual hmat_t<dcmplx> Overlap(const robs_t<dcmplx>* bs) const override
    {
        auto pw=dynamic_cast<const BasisSet::Band_FT_IBS*>(bs);
        assert(pw && "OrthoScalarFitter::Overlap requires a Band_FT_IBS (plane-wave) basis");
        return pw->Overlap(itsMap);
    }

    virtual void ReScale(double factor) override {for (auto& kv : itsMap) kv.second *= factor;}

    //! ScalarFunction (core): the fitted POTENTIAL v_xc,fit(r) = Re Σ_dm V-tilde(dm) e^{i(B·dm)·r}, via the
    //! basis's G_FieldEvaluator (DIP) -- what the GUI plots; the seed of the v_xc - v_xc,fit fit-residual.
    virtual double  operator()(const rvec3_t& r) const override {return FieldEval().EvalField(itsMap, r);}
    virtual rvec3_t Gradient  (const rvec3_t& r) const override {return FieldEval().EvalFieldGradient(itsMap, r);}

    virtual std::ostream& Write(std::ostream& os) const override
        {return os << "OrthoScalarFitter (orthonormal G-space overlap fit)" << std::endl;}

private:
    //! Reach the fit basis's inverse-transform capability (the DIP seam) from the neutral face we hold.
    const BasisSet::G_FieldEvaluator& FieldEval() const
    {
        auto* fe=dynamic_cast<const BasisSet::G_FieldEvaluator*>(itsFitBasis.get());
        assert(fe && "OrthoScalarFitter: the {G} fit basis must provide G_FieldEvaluator to evaluate v_xc,fit(r)");
        return *fe;
    }
    fbs_t      itsFitBasis;   //!< the {G} fit basis (the factory seam; inert until the denser-{G} upgrade)
    ΔG_Map     itsMap;        //!< the fit = the potential's V-tilde (received in DoFit)
};

} //namespace
