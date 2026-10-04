//! \file GPWTolerances.C
//! \brief The numerical tolerances of the Gaussian-plane-wave basis, as ONE typed value (D-ENV step 5).
//!
//! These were environment variables read inside the evaluator (a different process-wide static each).  They are tier-2
//! numerical-method policy: sensible defaults, almost never touched, but a CP2K-parity run or a convergence study sets them,
//! so they must be reproducible from a run description and not from a shell.  The facade carries one of these in its options
//! and hands it to the evaluator at construction.  The environment is NOT a way to set them (D-ENV step 6a): the input deck is
//! (`solid.tolerances.*`, or `--set`); the old GPW_* variables are retired, ignored and reported (qchem::WarnRetiredEnvironment).
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
    //! measured sub-mHa against CP2K.
    double vlocEps = 1.0e-5;
    //! The ABSOLUTE pair->level rule kappa (Ha per unit pair exponent) of the local-PP sweeps: every pair's spectral tail is
    //! bounded by e^{-kappa/2} independent of the field's sharpness.  30 Ha (e^{-15}) is CP2K's REL_CUTOFF 60 Ry; 60 is the
    //! self-convergence check.
    double localPPRelCutoff = 30.0;
    //! The HartreeOnly raster routing floor, as a FRACTION of alpha_max (only used with \c RasterFields::HartreeOnly).  beta=0
    //! (pure pair bandwidth) DIVERGES (+904 Ha), so the floor protects the density too.
    double relFieldSharp = 1.0/3.0;
    //! EXPLICIT multigrid cutoff list (Ha, descending) replacing the factor-ladder AND the top completion rung: level 0 stays
    //! the reference grid, these follow.  For matching CP2K's `CUTOFF/3^i` ladder, which the factor-4 default cannot
    //! reproduce.  Empty = the automatic ladder.   Deck: `"mgridEcuts":[53.33,17.78,5.926]`.
    std::vector<double> mgridEcuts = {};
    //! Magnitude screen of the ANALYTIC 1E lattice sums (S, <p^2>, V_local) and of the collocation boxes' reach: a pair whose
    //! Gaussian product is below this is dropped (CP2K's EPS_PGF_ORB analogue).  Drops only sub-eps terms, so S stays PSD.  Applied to the pair-loop evaluator by \c Periodic_Gaussian_IBS::ApplyTolerances.
    double screenEps = 1.0e-10;
    //! KS-field core exponent as a FRACTION of alpha_max (the raster sharpness; pin by rho_lost/N, not wall-clock).
    double fieldSharp = 2.0/3.0;
    //! Forces the ABSOLUTE pair->level rule at this kappa (Ha) when > 0 (CP2K REL_CUTOFF matching); 0 = the evaluator's own
    //! rule.
    double relCutoff = 0.0;
    //! The COLLOCATION tolerance floor: sizes the (shell pair, offset) task list and is the lowest any screener may answer.
    //! DECOUPLED from \c screenEps -- a collocated density and an analytic lattice sum are not converged by one tolerance.
    double densityEps = 1.0e-10;

    //! One line naming every value that differs from the default (empty when none) -- what a run banner prints.
    std::string Describe() const;
    bool operator==(const GPWTolerances&) const = default;
};

} // namespace qchem::BasisSet::Gaussian
