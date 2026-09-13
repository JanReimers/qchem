// File: DensityMixer.C  The DENSITY-FACE of SCF convergence (doc/SCFStrategyPlan.md).
//
// A density mixer folds a freshly diagonalised density (rho_out) into the running one (rho_in) and returns
// the convergence gate ‖Δρ‖.  It is one of the four SCF role-seams: the SCFIterator owns the density
// LIFECYCLE (SetWorkingCD / lineage / GetTotalEnergy) and the mixer owns the POLICY + arithmetic + state.
// This module is the whole CLIENT SURFACE -- the abstract face and the factories.  The concrete mixers live
// in their own modules under Internal/ (linear D-mixing; Kerker and Pulay on the G-space field; the polarized
// channel composition) and nothing outside this library names them: a client asks a factory for a mixer and
// talks to it through tDensityMixer alone.
module;
#include <memory>
export module qchem.ChargeDensity.DensityMixer;
export import qchem.ChargeDensity;                 // tChargeDensity<T>, tDM_CD<T>, tMixableDensity<T>
import qchem.BasisSet;                             // tBasisSet<T> (the periodic factory's fit-basis source)
import qchem.Structure;                            // Structure (the periodic factory's cell)
import qchem.Types;                                // dcmplx

export namespace qchem::ChargeDensity
{

//! Signals available to an adaptive density-mixing policy once the energy + orbital gradient are known.
struct MixSignals { double E=0.0, FD=0.0, FDold=0.0; };

//! The density-face of SCF convergence.  The SCFIterator calls, per fixed-point iteration:
//!   working = fresh diagonalised density (already made the lineage head by the iterator)
//!   dρ = Mix(working, old);            // fold rho_out into the running density
//!   ... iterator computes E, [F,D] ...
//!   if (WantsReDamp(sig)) { iterator reseats working->fresh; dρ = ReDampMix(working, old); recompute E; }
//!   UpdateRelax(sig);
//! and drives the next Fock from FockDensity(working).  A non-adaptive mixer takes the no-op defaults for the
//! three adaptive hooks.
template <class T> class tDensityMixer
{
public:
    //! THE SUBJECT OF MIXING: a density that can be mixed in place -- NOT the whole tDM_CD contract (ISP;
    //! see tMixableDensity).  A plain reference, not a shared_ptr: mixing never re-seats the caller's
    //! pointer (it mutates the pointee or keeps its own running ρ̃), and the caller necessarily outlives
    //! the call, so ownership never needed to be shared with the mixer at all.
    typedef tMixableDensity<T> cd_t;
    virtual ~tDensityMixer() {}
    //! Fold the fresh \a working density into the running one; returns ‖Δρ‖ (the convergence gate).
    virtual double Mix(cd_t& working, const cd_t& old) = 0;
    //! The density that drives the NEXT Fock: the working density itself for a D-mixer, the running mixed
    //! field for a G-space mixer.
    virtual const tChargeDensity<T>* FockDensity(const cd_t& working) const { return &working; }
    //! \brief Deposit the DM-backed density the next mixed field is built FROM, so a quadrature consumer can
    //! reach the EXACT density through \c tDM_Sourced_CD while Hartree keeps the preconditioned field.
    //!
    //! SHARED, and a SEPARATE hook rather than a widened \c Mix: \c Mix's subject is deliberately a plain
    //! reference (mixing never re-seats the caller's pointer and never outlives the call), and that reasoning
    //! is still right for mixing.  What is retained HERE outlives the call by construction -- XC samples it
    //! later in the same iteration -- so it is a different question and gets its own signature rather than
    //! quietly changing what \c Mix's argument means.
    //!
    //! No-op by default, which is correct for every D-mixing mixer: its \c FockDensity already IS the
    //! DM-backed density, so there is nothing to reach around.
    virtual void SetDMSource(std::shared_ptr<const tDM_CD<T>>) {}
    //! The current step size α (for the SCF trace only).
    virtual double GetRelax() const = 0;
    //! \brief The step size ACTUALLY DELIVERED, \f$\alpha_{\rm eff}\f$ -- for a preconditioned mixer the
    //! fraction of the update that survived the filter, which is what the SCF is really stepping at.
    //! Defaults to \c GetRelax(), exactly right for an unpreconditioned (linear) mix where the two coincide.
    //! Reported beside α in the ρ_mix column so a user can SEE the preconditioner working rather than infer
    //! it: on MnO it falls from 0.33 to 0.20 as the residual migrates into the damped low-G band, which is
    //! the difference between a healthy run and a stalling one and was previously invisible.
    virtual double EffectiveRelax() const { return GetRelax(); }
    //! A 3-char self-identifier for the per-iteration ρ_mix column (doc/GPWPlan1.md item 2).
    virtual const char* Tag() const { return "Lin"; }
    //! Adaptive [F,D]-keyed policy (the D-mixer's; no-op elsewhere).  Post-energy re-damp on divergence.
    virtual bool   WantsReDamp(const MixSignals&) const { return false; }
    //! Re-mix \a working (already reseated to the fresh density by the iterator) more aggressively; ‖Δρ‖.
    virtual double ReDampMix(cd_t& /*working*/, const cd_t& /*old*/) { return 0.0; }
    //! Grow/clamp the step for the next iteration.
    virtual void   UpdateRelax(const MixSignals&) {}
};

//=========================================================================================================
//  FACTORIES.  The implementations (all in Internal/, none of them a client concern):
//    * LinearMixer<T>          -- ρ_next = (1−α)ρ_in + α ρ_out on the density MATRIX (IrrepCD::MixIn), plus the
//                                 [F,D]-keyed adaptive α.  α=1 is passthrough, so there is NO NullMixer: "no
//                                 mixing" == the molecular default (StartingRelaxRo defaults to 1.0).
//    * KerkerMixer             -- Kerker-preconditioned ρ̃(G) mixing on the G-space field (periodic / dcmplx).
//    * PulayMixer              -- density-DIIS over a ρ̃ history, Kerker-preconditioned (periodic / dcmplx).
//    * PolarizedDensityMixer   -- one G-space leaf per channel ((ρ↑,ρ↓) or (ρ,m)) with ONE joint history.
//=========================================================================================================

//! The structure-neutral density mixer: plain linear D-mixing.  \a relax0 = StartingRelaxRo
//! (α=1 => passthrough).  Any run can use this one -- it asks nothing of the geometry.
template <class T> std::unique_ptr<tDensityMixer<T>> MakeLinearMixer(double relax0);

//! \brief The PERIODIC G-space mixer: Pulay when \a pulayDepth>0, else Kerker; one per SPIN CHANNEL on a
//! polarized run (PolarizedDensityMixer).
//!
//! Only a solid run asks for this, and a solid run HAS the periodic pieces by construction -- so the
//! Orbital_DFT_IBS<dcmplx> basis / UnitCell / FourierDensity faces are PRECONDITIONS here, not things to probe for.
//! This used to be one \c MakeDensityMixer that ran a three-way capability probe and fell back to linear
//! D-mixing with a warning: a periodic-vs-molecular decision sitting one layer too low.  The caller that
//! KNOWS (\c SolidSCFIterator::CreateMixer) now makes it, and a violated precondition THROWS rather than
//! silently mixing the wrong way for a whole run.  See doc/CleanupCandidates.md V1.10b.
//!
//! NB the SPIN branching below is a different question and keeps its graceful fallback: a spin-resolved
//! density that cannot hand out mutable channels is a real configuration to degrade from, not a broken
//! precondition -- so it warns and takes linear D-mixing, which at least keeps both channels.
//! \a cuspDeficit (N4, default false) picks the mixer that ALSO deposits the cusp-deficit correction for
//! \f$V_{xc}\f$.  It is a factory decision on purpose: CP2K parity is a property of WHICH MIXER WAS BUILT,
//! so the plain Kerker mixer stays bit-identical and no consumer has to read a flag to get it.
std::unique_ptr<tDensityMixer<dcmplx>> MakePeriodicMixer(
    double relax0, double kerkerG0, int pulayDepth, int pulayStart,
    const BasisSet::tBasisSet<dcmplx>* basis, const Structure* structure, const tDM_CD<dcmplx>* seed,
    bool cuspDeficit=false);

} //namespace
