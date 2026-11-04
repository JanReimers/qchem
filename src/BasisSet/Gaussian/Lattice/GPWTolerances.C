//! \file GPWTolerances.C
//! \brief The numerical tolerances of the Gaussian-plane-wave basis, as ONE typed value (D-ENV step 5).
//!
//! These were environment variables read inside the evaluator (a different process-wide static each).  They are tier-2
//! numerical-method policy: sensible defaults, almost never touched, but a CP2K-parity run or a convergence study sets them,
//! so they must be reproducible from a run description and not from a shell.  The facade carries one of these in its options
//! and hands it to the evaluator at construction.  The environment is NOT a way to set them (D-ENV step 6a): the input deck is
//! (`solid.tolerances.*`, or `--set`); the old GPW_* variables are retired, ignored and reported (qchem::WarnRetiredEnvironment).
module;
#include <optional>
#include <string>
#include <vector>
export module qchem.BasisSet.Gaussian.Lattice.GPWTolerances;

export namespace qchem::BasisSet::Gaussian
{

//! \brief The V_loc long routing tolerance a run WITHOUT an accuracy target uses (\c GPWTolerances::vlocEps 0): e^{-11.5}-class tails,
//! measured sub-mHa against CP2K -- the value every run used before it was derived from the target.
inline constexpr double kVlocEpsNoTarget = 1.0e-5;

//! \brief The GPW evaluator's tolerances.  Every default is today's behaviour; none of them is a convenience knob.
struct GPWTolerances
{
    //! Tolerance of the local-pseudopotential LONG-range G-ball (the harmonic per-pair rule that routes each pair over the multigrid
    //! levels).  0 = DERIVED (2026-11-02): the facade writes \c GridAccuracy::VlocEps() (the rigorous per-pair map corrected by the
    //! measured \c GridAccuracy::vlocRoutingEta, expert) when the run has an accuracy target, and the evaluator falls back to \c kVlocEpsNoTarget without one (stated rasters,
    //! accuracy 0).  A value > 0 given by the deck/caller is STATED: it WINS over the derived one and is honoured as is.
    double vlocEps = 0.0;
    //! REPORT-ONLY, written with a DERIVED \c vlocEps (0 = stated, or the no-target fallback): its two pieces, the rigorous per-pair map
    //! (\c GridAccuracy::RigorousVlocEps) and the measured calibration it was divided by (\c GridAccuracy::vlocRoutingEta).  The banner and
    //! the grids report print rigorous, correction and result; no number depends on them.
    double vlocEpsRigorous = 0.0;
    double vlocRoutingEta  = 0.0;
    bool VlocEpsDerived() const {return vlocRoutingEta>0.0;}
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
    //! \brief The RASTER demands (doc/GridEpsilonPlan.md §1 rows 1-2, C2 step 3), each an alias exponent x = G_N^2/p: every level's
    //! raster reaches G_N^2 >= x p_cap, p_cap = the tightest effective pair exponent the pair->level rule can route to that level, so
    //! the routing and the raster are ONE condition (a pair the rule sends to a level is sampled at x or better).  0 = the ball's own
    //! raster (2m+1: CP2K's arrangement, today's BallOnly).  They are DERIVED from the accuracy target by the facade
    //! (\c GridAccuracy: rigorous, stated calibration, or the run's own pilot audit) -- set one here only to state a raster.
    double rasterX     = 0.0;   //!< the density ladder (collocation + integrate-back: Hartree, and the uniform XC's input)
    double vlocRasterX = 0.0;   //!< V_loc long's own ladder (its static integrate-back; the term forces see first, M2c verdict 3)
    double xcRasterX   = 0.0;   //!< the reference level L0 again, for the UNIFORM XC raster (v_xc[rho] sampled there; MEASURED, no bound)

    //! One line naming every value that differs from the default (empty when none) -- what a run banner prints.
    std::string Describe() const;
    bool operator==(const GPWTolerances&) const = default;
};

//! \brief The alias exponent x = G_N^2/p each raster DELIVERS at the tightest pair it serves (C2 step 3) -- what the facade's pilot
//! audit consumes to fit a term's measured tail at the run's TRUE raster (the demand is a floor; the ball and the 5-smooth menu can
//! exceed it).  0 = that raster is not built (V_long before its sweep, or a model with no closed Gaussian long part).
struct RasterExponents
{
    double density  = 0.0;   //!< the finest density level at p = 2 alpha_max (Hartree's binding pair)
    double vlocLong = 0.0;   //!< V_loc long's finest level at p = 2 alpha_max
    double xcL0     = 0.0;   //!< the reference level L0 at the tightest pair routed to it (the uniform XC samples v_xc[rho] there)
};

//! \brief One raster's GRID-AUDIT verdict (the `grids.audit*` row's status and alias): what the facade's fixed-D grid audit records for
//! each side and bounds a raster term's delivered error at (the coarse side's alias).  Asked of the raster's owner
//! (\c GpwRasterLadder::DeliveredRasterVerdicts), never read back out of the run report (V-REPORT-WRITEONLY).
struct RasterVerdict
{
    std::string audit;          //!< the row's key: "auditDensity" | "auditVlocLong"
    std::string status;         //!< "OK" | "UNDER" (the ball tail above the tolerance it was built to)
    double      alias = 0.0;    //!< e^{-x}, x = G_N^2/p at the tightest pair the raster serves
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
//! measured by `gridAudit` (doc/GridEpsilonPlan.md §5 M2c): eta = delivered / bound <= 1, so a term meets eps when eta*Bound(x) <= eps,
//! i.e. the sizing uses L(eps/eta) -- LOOSER than rigorous, never tighter.  1 = rigorous mode (geometry/force
//! validation asks for it).  M2c found eta NOT flat across materials (V_long 1e-4 ... 0.15, Hartree 5e-6 ... 1e-2), so it is to be
//! SELF-CALIBRATED per run by a pilot audit, never a shipped per-class table (\c GridAccuracy::eta empty, the default).
//! For FORCES eta is judged on max|dH_ij|, never on energies (dE cancels inside Tr D dH: MnO V_long 0.26 Ha vs dE 5 uHa).
struct AccuracyCalibration
{
    double vlocLong = 1.0;   //!< eta of the V_loc long term
    double density  = 1.0;   //!< eta of the density-raster terms (Hartree; NOT the XC raster -- see GridAccuracy)
    bool operator==(const AccuracyCalibration&) const = default;
};

//! \brief WHAT THE RUN'S CLIENT CONSUMES of the Hamiltonian -- the measure the ONE accuracy target is judged on (user ruling
//! 2026-10-31, doc/GridEpsilonPlan.md §1): a single-point ENERGY is judged on the energy error (eps per atom), and anything that
//! produces forces, a geometry or a surface on the matrix elements (eps per element, max|dH_ij|: the energy error cancels inside
//! Tr D dH, M2c verdict 3).  The run ANSWERS it (who built it), the user never picks a metric.
enum class ClientConsumes { Energy, MatrixElements };

//! \brief The ONE accuracy target (deck `solid.accuracy`) and its calibration (`solid.accuracyCalibration`, expert).  C2 step 3: the
//! raster demands (\c GPWTolerances::rasterX ...) are DERIVED from it -- by the stated calibration when there is one (eta = 1 is the
//! rigorous mode), else by the run's own PILOT audit (self-calibrated: M2c found eta not flat across materials).
struct GridAccuracy
{
    double eps = 1.0e-5;   //!< Ha per atom on an energy run, Ha per matrix element on a geometry run (\c ClientConsumes)
    //! The STATED calibration (expert); none = SELF-CALIBRATED by a pilot audit (the default).
    std::optional<AccuracyCalibration> eta;
    //! L_rigorous (eta = 1) and L_calibrated = L(eps/eta) (eta = delivered/bound <= 1 loosens) for one tail.  Clamped at x = n (the
    //! tail's maximum) when eps/eta exceeds the bound's maximum: no raster beyond the floor is needed for that term.
    double Rigorous(const AliasTail& t) const {return t.Exponent(eps);}
    double Calibrated(const AliasTail& t, double etaTerm) const;
    //! The target one term's delivered error is held to: eps x the atom count on an energy run, eps per element on a geometry run.
    double Target(ClientConsumes c, size_t nAtoms) const {return c==ClientConsumes::Energy ? eps*double(nAtoms) : eps;}
    //! \brief RIGOROUS: the per-pair RELATIVE bound.  MakeLocalPPLong's harmonic rule (lnE = -ln vlocEps, ball edge 2 p_h lnE) routes each
    //! (pair, offset) to a level whose spectral tail, RELATIVE to that pair's own contribution, is bounded by vlocEps -- no S prefactor,
    //! no x^3 Hermite factor (S_Vlong = 2.35 is the RASTER alias tail's prefactor, a different error).  Turning that into Ha needs two
    //! things the bound does not give: a TERM SCALE (the pair's matrix-element magnitude) and the PAIR COUNT (an element, and a cell
    //! energy, sum many pairs and offsets whose truncation errors may add or cancel).  Neither is provable, so the rigorous map is the
    //! identity: the per-pair relative bound = eps.  \c vlocRoutingEta carries the scale and the pair sum together, MEASURED.
    double RigorousVlocEps() const {return eps;}
    //! \brief EMPIRICAL: the V_loc long routing calibration eta = delivered / bound (expert `solid.accuracyCalibration.vlocRouting`); it
    //! carries the term scale and the pair sum the rigorous piece cannot.  The default 10 is NOT itself measured; it decomposes as
    //! (measured max delivered/bound) x (margin), rounded up:
    //!   - measured 2026-11-02 against a reference that re-routes every pair (the pilot's ball audit), at vlocEps 1e-5 and 1e-7 alike:
    //!     max 6 per cell energy (MnO, 4 atoms; Si 2 atoms 1.4) and max 1.7 per matrix element (Si Gamma, geometry; low 1.0);
    //!   - times a >= 1.5 x margin (6 x 1.5 = 9), rounded to 10.
    //! A re-measurement updates the first factor; the margin is the policy.  The pilot's honest gate (every audited term against the
    //! target) catches a cell where the pair sum is worse.  Without it (fixed vlocEps 1e-5) MnO's V_loc long delivered 6e-5 against a
    //! 4e-5 target, alternating in sign as the ball doubled (the multigrid levels re-place the pairs).  A STATED absolute
    //! \c GPWTolerances::vlocEps (> 0) wins over both pieces.
    double vlocRoutingEta = 10.0;
    //! The routing tolerance the facade derives: the rigorous per-pair relative bound corrected by the measured calibration -- the ONE
    //! place the two pieces meet.
    double VlocEps() const {return RigorousVlocEps()/vlocRoutingEta;}
    //! \brief The STATED calibration's raster demands written into \a t (rasterX = the density tail's L(eps/eta), vlocRasterX = V_long's,
    //! xcRasterX = rasterX: the XC raster has no bound of its own).  THROWS when there is no stated calibration (that run self-calibrates).
    GPWTolerances Stated(GPWTolerances t) const;
    //! \brief HEURISTIC, the self-calibration's one step: a term that DELIVERED \a delivered at alias exponent \a xRun fixes its own scale
    //! C = delivered / (xRun^3 e^{-xRun}) under the rigorous SHAPE (\c AliasTail n = 3), and demands the x where C x^3 e^{-x} = \a target.
    //! 0 when it already meets the target (no demand beyond the run's raster).  eta = C / S_term is what the report prints.
    static double MeasuredDemand(double delivered, double xRun, double target);
    //! \brief MEASURED, for the uniform XC raster (no bound, no shape: v_xc[rho] is nonlinear, and Si's max|dH_xc| fell only as a POWER
    //! of x -- 1.0e-2 / 6.3e-3 / 1.1e-3 / 4.8e-4 at x 21 / 30 / 37 / 54, where the alias shape predicted e^{-x}): from two of the pilot's
    //! own rounds (x1, d1), (x2, d2), d ~ x^{-k} with k = ln(d1/d2)/ln(x2/x1), demand x2 (d2/target)^{1/k}.  0 when (x2, d2) meets the
    //! target, or when no decrease was measured (no slope: the caller falls back to \c MeasuredDemand).
    static double MeasuredPowerDemand(double x1, double d1, double x2, double d2, double target);
    //! \brief The run banner's accuracy line (after "[label run] accuracy: "): eps, the calibration (stated or self), per term class the
    //! rigorous x, and the raster demands USED.  The XC raster has NO rigorous bound (v_xc[rho] is nonlinear: it aliases even where
    //! rho~ is band-limited -- Si k222 434x over the rho-pair bound, M2c), so it is labelled MEASURED, never rigorous.
    std::string Describe(const GPWTolerances& used) const;
    bool operator==(const GridAccuracy&) const = default;
};

} // namespace qchem::BasisSet::Gaussian
