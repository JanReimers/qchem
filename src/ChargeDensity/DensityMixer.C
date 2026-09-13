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
//! and drives the next Fock from FockDensity(working).  A mixer whose step size ADAPTS to the loop's
//! signals additionally implements tAdaptiveMixer (below) -- a capability face, reached by cross-cast, so a
//! fixed-step mixer carries nothing for it.
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
    virtual const tChargeDensity<T>* FockDensity(const cd_t& working) const = 0;
    // (No SetDMSource here -- V1.18: which D a mixed field was built from is PROVENANCE the loop driver seats
    //  on the Fock density itself, through tDM_SourceSink.  A mixer neither stashes nor replays it.)
    //! The current step size α (for the SCF trace only).
    //! (No "effective α" beside it any more -- user ruling 2026-09-13: the fraction of a step that survived a
    //! preconditioner is not physics, nothing consumes it, and printing it implied something did.)
    virtual double GetRelax() const = 0;
    //! A 3-char self-identifier for the per-iteration ρ_mix column (doc/GPWPlan1.md item 2).
    virtual const char* Tag() const = 0;
};

//! \brief Capability face: a mixer whose STEP SIZE adapts to the loop's signals.  Only the linear D-mixer has
//! it (its [F,D]-keyed α), so it lives here and not on tDensityMixer -- the iterator cross-casts ONCE when
//! it builds the mixer, and a fixed-step mixer implements nothing.
//!
//! This replaced three defaulted hooks (WantsReDamp / ReDampMix / UpdateRelax, V1.18) whose choreography
//! had the ITERATOR rebuilding ρ_out from the wave function so the mixer could re-mix it -- unnecessary,
//! because the mix is LINEAR: the re-damped density is reachable from the already-mixed one (see
//! LinearMixer::Adapt).  The energy recompute stays the iterator's: a mixer is not an energy service, it
//! just says whether it changed the density under the caller's feet.
template <class T> class tAdaptiveMixer
{
public:
    typedef tMixableDensity<T> cd_t;
    virtual ~tAdaptiveMixer() {}
    //! POST-ENERGY hook, once this step's E and [F,D] are known.  May RE-MIX \a working with a smaller step;
    //! returns true when it did, with \a dRho the new convergence gate -- the caller must then recompute
    //! whatever it derived from \a working.  Either way, sets the step for the NEXT iteration.
    virtual bool Adapt(const MixSignals&, cd_t& working, const cd_t& old, double& dRho) = 0;
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
