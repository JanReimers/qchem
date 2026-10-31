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
    //! Gaussian product is below this is dropped (CP2K's EPS_PGF_ORB analogue).  Drops only sub-eps terms, so S stays PSD.  Bound on the molecular basis by its builder (\c ToleranceBinding, D-ENV-TOL option A).
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

//! \brief RIGOROUS (doc/GridEpsilonPlan.md §1, ruled 2026-10-31: "the rigorous math aspect well isolated from the heuristic part"):
//! one term class's raster ALIAS-TAIL model, err(x) <= S x^n e^{-x} at x = G_N^2/p, fitted ONCE against the closed-form oracles
//! (`GPW.RasterAliasCalibration`, M2b: an s+p+d pair's Hermite prefactor makes n = 3) and gated by them.  Pure math, no fudge inside:
//! at eta = 1 it guarantees the single-pair worst case.  Valid for x >= 4 (below it the image is not a small tail).
struct AliasTail
{
    double S;   //!< the term's scale (Ha for a matrix-element bound; relative for the density)
    double n;   //!< the Hermite prefactor's power
    //! The bound at a raster's x.
    double Bound(double x) const;
    //! \brief The x a raster must reach for the bound to meet \a eps: the root L > n of S L^n e^{-L} = eps (the DECREASING branch --
    //! the bare ln(S/eps) under-sizes by ~7 units of x at 1e-5).  THROWS for eps >= the bound's maximum (no raster needed).
    double Exponent(double eps) const;
};
//! M2b's fitted tails (n = 3 is the gate; GPW_UT `RasterAliasCalibration` asserts err <= S x^3 e^{-x} with S 2.5 / 0.4).
inline constexpr AliasTail kVlocLongTail{2.35, 3.0};   //!< V_loc long <chi|V|chi> (Ha)
inline constexpr AliasTail kDensityTail {0.35, 3.0};   //!< rho~ at low G (relative): Hartree's and the raster XC's input

//! \brief HEURISTIC (the same ruling): how far the DELIVERED error undershoots the rigorous bound, ONE factor per term class,
//! measured by `gridAudit` (doc/GridEpsilonPlan.md §5 M2c).  The sizing uses L(eta*eps).  1 = rigorous mode (geometry/force
//! validation asks for it).  M2c found eta NOT flat across materials (V_long 1e-4 ... 0.15, Hartree 5e-6 ... 1e-2), so it is to be
//! SELF-CALIBRATED per run by a pilot audit, never a shipped per-class table; until that lands the default is the rigorous 1.
//! For FORCES eta is judged on max|dH_ij|, never on energies (dE cancels inside Tr D dH: MnO V_long 0.26 Ha vs dE 5 uHa).
struct AccuracyCalibration
{
    double vlocLong = 1.0;   //!< eta of the V_loc long term
    double density  = 1.0;   //!< eta of the density-raster terms (Hartree; NOT the XC raster -- see GridAccuracy)
    bool operator==(const AccuracyCalibration&) const = default;
};

//! \brief The ONE accuracy target (deck `solid.accuracy`, Ha per Hamiltonian matrix element) and its calibration
//! (`solid.accuracyCalibration`, expert).  C2 step 2: carried, derived and REPORTED; the sizing still uses today's knobs (step 3
//! makes them derived -- the anchor-moving step), so the [accuracy] line prints both what the rule demands and what is used.
struct GridAccuracy
{
    double              eps = 1.0e-5;   //!< today's tests' value (kappa = 30 Ha <=> L = 15 on the local-PP sweeps)
    AccuracyCalibration eta;
    //! L_rigorous (eta = 1) and L_used-by-the-rule (eta applied) for one tail.
    double Rigorous(const AliasTail& t) const {return t.Exponent(eps);}
    double Calibrated(const AliasTail& t, double etaTerm) const {return t.Exponent(etaTerm*eps);}
    //! \brief The run banner's accuracy line (after "[label run] accuracy: "): eps, eta, and per term class the raster x the rule demands, rigorous and calibrated, beside
    //! the knob values the sizing uses today.  The XC raster has NO rigorous bound (v_xc[rho] is nonlinear: it aliases even where
    //! rho~ is band-limited -- Si k222 434x over the rho-pair bound, M2c), so it is labelled MEASURED, never rigorous.
    std::string Describe(const GPWTolerances& used) const;
    bool operator==(const GridAccuracy&) const = default;
};

} // namespace qchem::BasisSet::Gaussian
