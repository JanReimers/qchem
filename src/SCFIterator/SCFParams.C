// File: SCFParams.C Parameters used for controling SCF iteration and convergence.
module;
#include <stdexcept>
export module qchem.SCFParams;
import qchem.Types;

namespace qchem {

// Defaults are the proven molecular-HF recipe (matches M_HF_U / the bridge): they converge a
// neutral molecule (e.g. water/dzvp) out of the box, so a facade caller can write
// `calc.Converge({.NMaxIter=60})` and override only what they care about.  Positional aggregate
// init (`{20,1e-4,...}`) and per-field assignment (`par.MinΔρPerElectron=...`) both still work unchanged --
// default member initializers do not disqualify the aggregate.
export struct SCFParams
{
    size_t NMaxIter        = 20;     //Max allowed number of iterations
    //! THE density convergence gate (D-MINDRHO, 2026-10-22; was \c MinΔρPerElectron, whose meaning depended on the mixer): converged
    //! when the step would still move fewer than this many electrons PER ELECTRON, \f$\int|\rho_{out}-\rho_{in}|\,d^3r/N\f$
    //! on the run's own real-space quadrature (per spin channel for a polarized run; see \c Δρmeasure).
    double MinΔρPerElectron = 1e-4;
    //! WHAT THE GATE MEASURES.  L1PerElectron (default) = the one intensive residual above, computed once in the SCF
    //  iterator, the same number for every mixer, representation and system (a field mixer's in/out fields are presented
    //  alike).  It REPLACED "MixerResidual" -- whatever the mixer reported: max|ρ̃_out(G)−ρ̃_in(G)| on the Kerker/Pulay
    //  path, Σ‖ΔD‖_F/N on the linear-D path -- which the calibration (2026-10-22) found ~faithful for linear mixing but
    //  50-120x LOOSE for Kerker/Pulay (runs "converged" with ~2 % of the electrons still moving).  MaxΔD = CP2K's EPS_SCF
    //  measure: max_ij|D_ij−D'_ij| between successive OUTPUT density matrices, un-normalised (qs_scf_loop_utils.F); then
    //  \c MinΔρPerElectron is compared with THAT -- select it with the EPS_SCF value when an iteration count is to be
    //  compared with a CP2K deck's.
    enum class Measure { L1PerElectron, MaxΔD };
    Measure Δρmeasure      = Measure::L1PerElectron;
    double MinΔFD          = 1e-7;   //Minimum delta in [F,D] (|[F,D]-[F,D]_old|) for convergence.
    double MinΔE           = 1e30;   //Minimum RELATIVE total-energy change |ΔE/E| for convergence (default
                                     //  off).  The physical gate for a NON-variational SCF (density-fit /
                                     //  GPW collocation) whose [F,D]/Δρ limit-cycle above the fit floor while
                                     //  the total energy is settled -- see doc/GPWPlan.md.
    double MinVirial       = 1e30;   //Max |error| in the virial ratio -V/K.  DEFAULT OFF (1e30), because
                                     //  the test below is |err| < MinVirial: a SMALL threshold is STRICT,
                                     //  not lax.  The old 1e-13 default read as "off" but was the tightest
                                     //  possible gate -- unsatisfiable, so it could only ever BLOCK
                                     //  convergence.  It never bit because every serious caller sets this.
                                     //  Turn it ON (e.g. 4e-2) for what it is actually good at: very
                                     //  high-precision ATOM work on BSpline bases, where -V/K=2 is a real
                                     //  accuracy check (user 2026-08-10).  Meaningless under a
                                     //  pseudopotential -- V1.27's IsVirialValid() drops the gate there
                                     //  automatically, so no caller has to remember.
    double MinFD           = 1e-5;   //Minimum error from SCF accelerator.  i.e. [F,D] (orbital gradient).
    double StartingRelaxRo = 1.0;    //relaxation for mixing Ro.  Dynamically adjusted during iterations.
    double MergeTol        = 1e-4;   //Merge eigen levels (like Px,Py Pz) that are equal within +/- MergeTol
    bool   Verbose         = true;   //Display iteration details.  ON by default (user, 2026-10-05): a run's console without its iteration table hides where the time goes.
    //! \brief Feed \f$V_{xc}\f$ the CUSP-DEFICIT correction \f$\rho[D]-\tilde\rho_{fit}\f$ on top of the mixed density,
    //! instead of the fitted \f$\tilde\rho_{mix}\f$ alone (doc/OpenWork.md N4) -- how a Kerker/Pulay run meets doc/Pins.md pin 28
    //! (XC from ρ[D]): at the fixed point \f$\tilde\rho_{mix}=\tilde\rho_{fit}\f$ and XC sees ρ[D] exactly.
    //!
    //! A MIXER-CONSTRUCTION policy, deliberately, not a switch the XC code reads: CP2K parity is a property
    //! of WHICH MIXER WAS BUILT, so it must be decidable at the factory.  TRUE BY DEFAULT since 2026-10-22 (pin 28; user);
    //! false = the pre-pin-28 mixer, XC on the fitted ρ̃_mix (the correction is not even formed) -- a control, not a recipe.
    //! Read only by the Kerker/Pulay mixers (linear D-mixing and GDM already feed XC ρ[D]).  Wired on Kerker and Pulay,
    //! unpolarized and in both spin coordinates; the pair (raster) XC route serves it with its raw shadow and refuses it on
    //! the BALL fallback (DensitySampler_Pair.C).
    bool   XCCuspDeficit   = true;
    double KerkerG0        = 0.0;    //Kerker density-mixing screening wavevector G0 (a.u.^-1).  0 (default) =
                                     //  OFF -> classic linear density-matrix mixing (atoms/molecules unchanged).
                                     //  >0 -> preconditioned rho-mixing on the periodic (GPW/PW) path: damp the
                                     //  low-G charge-transfer slosh by G^2/(G^2+G0^2) (breaks the ionic-crystal
                                     //  limit cycle -- see doc/GPWPlan.md).  Typical G0 ~ 1.0.  dcmplx path only.
    bool   UseMOM          = false;  //Maximum Overlap Method: occupy the orbitals with the largest overlap onto a
                                     //  FIXED reference occupied subspace instead of the lowest eigenvalues.  OFF
                                     //  (default) -> plain aufbau (atoms/molecules unchanged).  ON -> occupied-
                                     //  subspace continuity, the fix for a diving diffuse virtual being aufbau-
                                     //  captured (the NaF Γ occupation-swap instability -- see doc/GPWPlan §0b″).
    int    PulayDepth      = 0;      //Density-face Pulay (density-DIIS) history depth on the periodic ρ̃ path.  0
                                     //  (default) = plain Kerker/linear mixing.  >0 = Kerker-preconditioned Pulay
                                     //  mixing keeping this many (ρ̃_in, ρ̃_out) pairs (VASP/QE/CP2K scheme; the
                                     //  quasi-Newton cure for the charge-transfer slosh -- see doc/SCFStrategyPlan.md).
    int    PulayStart      = 0;      //Prime with plain Kerker for this many iterations before engaging Pulay (only
                                     //  used when PulayDepth>0).  History-based mixing is UNSTABLE far from the fixed
                                     //  point (residuals not yet in the linear-response regime), so descend first,
                                     //  then accelerate -- the density-side of the ladder hand-off.  0 = immediate.
    int    MOMStartIter    = 10;     //Delayed-IMOM reference-capture iteration (only used when UseMOM).  Run plain
                                     //  aufbau for this many fills so the SCF descends to the physical fixed point,
                                     //  THEN capture the occupied subspace ONCE and hold it fixed (the seed is
                                     //  mid-transient, so capturing at iteration 0 anchors garbage).
    double SmearingkT      = 0.0;    //Fermi-Dirac smearing electronic temperature kT (Hartree).  0 (default) = OFF
                                     //  -> integer aufbau (atoms/molecules/gapped solids unchanged).  >0 -> fractional
                                     //  occupations f_i=1/(1+exp((ε_i−μ)/kT)) with μ solved per block by bisection on
                                     //  Σg_i f_i=nₑ; the Mermin −TS enters ONLY the free energy (EnergyBreakdown::MinusTS
                                     //  -> GetTotalEnergy becomes A=E−TS), never H.  The cure for near-gapless occupation
                                     //  flapping (NaF Ecut=160; doc/GPWPlan1.md 4b).  Typical kT ~ 1e-3..1e-2 Ha.
                                     //  Composes with UseMOM via MOMSmearPenalty below (MOM-masked Fermi).
    double MOMSmearPenalty = 0.0;    //Advanced (only when UseMOM AND SmearingkT>0): the MOM-overlap penalty Λ that
                                     //  pushes LOW-overlap states (diffuse ghosts) up in EFFECTIVE energy
                                     //  ε_i+Λ(1−s_i)² for the Fermi fill, so they stay empty by CHARACTER while
                                     //  the physical frontier smears by TRUE energy.  0 (default) = smearing alone
                                     //  (fills by energy; a diving diffuse ghost would be occupied).  Λ needs only
                                     //  exceed the ghost's dive depth (Ha-scale); doc/GPWPlan1.md 4b.

    bool StopOnAccelExhausted = false; //End the SCF loop when the accelerator declares itself Exhausted()
                                     //  (the LADDER's hand-off trigger fired but the next rung is NOT
                                     //  ENGAGEABLE -- e.g. DIIS out of steam and GDM vetoed by smeared
                                     //  occupations; doc/SCFStrategyPlan.md "Ladder exhaustion bail-out",
                                     //  run 49's 30 wasted DIIS iterations).  OFF (default) = keep grinding
                                     //  (run 38's post-veto tail DID eventually converge).  The annealed
                                     //  driver turns this ON for NON-FINAL stages: their job -- get the
                                     //  density near the next stage's fixed point -- is done at the tie
                                     //  floor, and the cold stage re-arms GDM.

    //! ADVANCED (doc/GPWPlan1.md item 1): the 0h occupation-guard tunings.  Sensible defaults that only
    //! rarely need touching -- the guard exists precisely so users do NOT hand-tune around a pinned
    //! excited state (the expert-system direction: detect, self-correct, and never be silent).
    struct MOMGuardParams
    {
        int HolePersistence = 3;     //Consecutive non-aufbau (hole) iterations before the MOM reference is
                                     //  released (drop + re-arm the delayed-IMOM capture).
        int MaxReleases     = 2;     //Releases per run; a hole surviving them is WARNED loudly at run end.
    };
    MOMGuardParams Guard;            //Advanced: the MOM guard (only consulted when UseMOM).

    //! \brief WHICH density mixer (OpenWork §3 "Kerker/Pulay as the DEFAULT mixer"; user 2026-10-04: *"in general
    //! Kerker/Pulay should be the default mixer"*).
    //!
    //! \c Auto (default) is resolved by the ITERATOR, which knows whether it can mix in G space and whether its accelerator
    //! keeps a Fock history (\c ResolveMixing below): a periodic run with NO Fock extrapolator gets Pulay on D (no Kerker:
    //! user 2026-10-26, pin 28) with ITS OWN recipe (\c theAutoMixing -- \c StartingRelaxRo's 1.0 default is LINEAR's step);
    //! a molecule/atom, or a run beside DIIS, gets linear D-mixing at \c StartingRelaxRo (Pulay beside DIIS fights; the solid
    //! facade's \c Auto accelerator is the other half of the pair).  A caller
    //! that already set a G-space knob (\c KerkerG0>0, \c PulayDepth>0) keeps exactly that recipe under \c Auto, so
    //! every hand-tuned recipe in the tree means what it meant.  The explicit modes say what they read: \c Linear =
    //! \c StartingRelaxRo only (the G-space knobs are ignored); \c Kerker = \c KerkerG0 (>0), no history; \c Pulay =
    //! \c PulayDepth (>0), \c PulayStart, \c KerkerG0 (0 = an unpreconditioned history).  Kerker and Pulay THROW on a
    //! run that cannot mix in G space.
    enum class Mixing { Auto, Linear, Kerker, Pulay };
    Mixing Mixer           = Mixing::Auto;
};

//! The density-mixing recipe \c Mixing::Auto resolves to on a periodic run (\c ResolveMixing).
export struct MixingRecipe
{
    double relax;   //!< the mixing step α (→ \c StartingRelaxRo)
    double G0;      //!< Kerker screening wavevector, a.u.⁻¹ (→ \c KerkerG0)
    int    depth;   //!< Pulay history depth (→ \c PulayDepth)
    int    start;   //!< plain-Kerker priming iterations before the history engages (→ \c PulayStart)
};
//! THE default periodic mixer: PULAY ON D (user 2026-10-26) -- the history on density matrices, so V_xc sees rho[D] with no
//! fit and no cusp deficit (pin 28), and no rho-tilde is built per iteration (the Kerker arm cost a vacuum box 20x per
//! iteration).  The step/history are the banked Kerker/Pulay recipe's (free MnO AFM-II, Si k211 FD gate), close to VASP's
//! AMIX 0.4 and CP2K's ALPHA 0.4 + NBUFFER; Kerker (G0 > 0) is a stated choice, not the default.
export inline constexpr MixingRecipe theAutoMixing{.relax=0.45, .G0=0.0, .depth=8, .start=5};

//! \brief Resolve \c p.Mixer to a concrete mode with its knobs filled; the result's \c Mixer is never \c Auto.
//! \a gSpace = can this run mix in G space (the periodic iterator asks its density); \a fockExtrapolates = the run's accelerator
//! keeps a Fock history (DIIS, a ladder with a DIIS rung), beside which \c Auto stays LINEAR -- a density history next to a Fock
//! history fights (measured 2026-10-22); \a recipe = what \c Auto becomes otherwise.
//! Idempotent: resolving a resolved value returns it unchanged.  Throws on an explicit G-space mode without its knob,
//! or without \a gSpace.
export inline SCFParams ResolveMixing(const SCFParams& p, bool gSpace, bool fockExtrapolates=false,
                                     const MixingRecipe& recipe=theAutoMixing)
{
    using M = SCFParams::Mixing;
    SCFParams r=p;
    switch (p.Mixer)
    {
    case M::Linear:
        return r;
    case M::Kerker:
    case M::Pulay:
        if (!gSpace && !(p.Mixer==M::Pulay && p.KerkerG0<=0.0))   // Pulay at G0 = 0 mixes D: any run (PulayDMixer)
            throw std::invalid_argument("SCFParams: Kerker / Kerker-preconditioned Pulay mixing needs a periodic (G-space) run; "
                                        "use Mixing::Linear, or Mixing::Pulay with KerkerG0 = 0 (Pulay on D).");
        if (p.Mixer==M::Kerker && p.KerkerG0<=0.0)
            throw std::invalid_argument("SCFParams: Mixing::Kerker needs KerkerG0 > 0.");
        if (p.Mixer==M::Pulay && p.PulayDepth<=0)
            throw std::invalid_argument("SCFParams: Mixing::Pulay needs PulayDepth > 0.");
        if (p.Mixer==M::Kerker) r.PulayDepth=0;
        return r;
    case M::Auto:
        break;
    }
    if (!gSpace)              { r.Mixer=M::Linear; return r; }   // molecules / atoms: today's linear D-mixing
    if (p.PulayDepth>0)       { r.Mixer=M::Pulay;  return r; }   // a hand-set recipe is kept as it was
    if (p.KerkerG0>0.0)       { r.Mixer=M::Kerker; return r; }
    if (fockExtrapolates)     { r.Mixer=M::Linear; return r; }   // beside DIIS: today's linear + DIIS
    r.StartingRelaxRo=recipe.relax; r.KerkerG0=recipe.G0;
    r.Mixer=M::Pulay; r.PulayDepth=recipe.depth; r.PulayStart=recipe.start;
    return r;
}



} // namespace qchem