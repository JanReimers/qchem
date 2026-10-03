//! \file GPWTolerances.C
//! \brief The numerical tolerances of the Gaussian-plane-wave basis, as ONE typed value (D-ENV step 5).
//!
//! These were environment variables read inside the evaluator (a different process-wide static each).  They are tier-2
//! numerical-method policy: sensible defaults, almost never touched, but a CP2K-parity run or a convergence study sets them,
//! so they must be reproducible from a run description and not from a shell.  The facade carries one of these in its options
//! and hands it to the evaluator at construction; the environment survives only as an OVERRIDE LAYER applied where the
//! options are resolved (\ref qchem::BasisSet::Gaussian::ApplyEnvOverrides), and what an override changed is reported.
module;
#include <string>
#include <vector>
export module qchem.BasisSet.Gaussian.Lattice.GPWTolerances;

export namespace qchem::BasisSet::Gaussian
{

//! \brief The GPW evaluator's tolerances.  Every default is today's behaviour; none of them is a convenience knob.
struct GPWTolerances
{
    //! Tolerance of the local-pseudopotential LONG-range G-ball (the harmonic per-pair rule): e^{-8.5}-class tails,
    //! measured sub-mHa against CP2K.  Env override `GPW_VLOC_EPS`.
    double vlocEps = 1.0e-5;
    //! The ABSOLUTE pair->level rule kappa (Ha per unit pair exponent) of the local-PP sweeps: every pair's spectral tail is
    //! bounded by e^{-kappa/2} independent of the field's sharpness.  30 Ha (e^{-15}) is CP2K's REL_CUTOFF 60 Ry; 60 is the
    //! self-convergence check.  Env override `GPW_LOCALPP_RELCUTOFF`.
    double localPPRelCutoff = 30.0;
    //! The HartreeOnly raster routing floor, as a FRACTION of alpha_max (only used with \c RasterFields::HartreeOnly).  beta=0
    //! (pure pair bandwidth) DIVERGES (+904 Ha), so the floor protects the density too.  Env override `GPW_RELFIELDSHARP`.
    double relFieldSharp = 1.0/3.0;
    //! EXPLICIT multigrid cutoff list (Ha, descending) replacing the factor-ladder AND the top completion rung: level 0 stays
    //! the reference grid, these follow.  For matching CP2K's `CUTOFF/3^i` ladder, which the factor-4 default cannot
    //! reproduce.  Empty = the automatic ladder.  Env override `GPW_MGRID_ECUTS="53.33,17.78,5.926"`.
    std::vector<double> mgridEcuts = {};
    //! Magnitude screen of the ANALYTIC 1E lattice sums (S, <p^2>, V_local) and of the collocation boxes' reach: a pair whose
    //! Gaussian product is below this is dropped (CP2K's EPS_PGF_ORB analogue).  Drops only sub-eps terms, so S stays PSD.
    //! Env override `GPW_SCREEN_EPS`.  Applied to the pair-loop evaluator by \c Periodic_Gaussian_IBS::ApplyTolerances.
    double screenEps = 1.0e-10;
    //! KS-field core exponent as a FRACTION of alpha_max (the raster sharpness; pin by rho_lost/N, not wall-clock).
    //! Env override `GPW_FIELDSHARP`.
    double fieldSharp = 2.0/3.0;
    //! Forces the ABSOLUTE pair->level rule at this kappa (Ha) when > 0 (CP2K REL_CUTOFF matching); 0 = the evaluator's own
    //! rule.  Env override `GPW_RELCUTOFF`.
    double relCutoff = 0.0;
    //! The COLLOCATION tolerance floor: sizes the (shell pair, offset) task list and is the lowest any screener may answer.
    //! DECOUPLED from \c screenEps -- a collocated density and an analytic lattice sum are not converged by one tolerance.
    //! Env override `GPW_DENSITY_EPS`.
    double densityEps = 1.0e-10;

    //! One line naming every value that differs from the default (empty when none) -- what a run banner prints.
    std::string Describe() const;
    bool operator==(const GPWTolerances&) const = default;
};

//! \brief Apply the environment as the OVERRIDE LAYER onto \a t (the environment wins; each override is named in the
//! returned string so the banner can say it).  Called ONCE, where the facade resolves its options -- never from a leaf.
std::string ApplyEnvOverrides(GPWTolerances& t);

} // namespace qchem::BasisSet::Gaussian
