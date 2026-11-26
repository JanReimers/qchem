// File: Calculation/RasterPilotXC.C  The raster pilot's demand rule for the UNIFORM XC term (lane A, 2026-11-25; OpenWork "the raster pilot
// never weighs a bigger ball against a finer raster, and its XC demand underestimates every jump").
//
// v_xc[rho] on the uniform route is sampled on ladder level L0, so the XC term's delivered error can come from two rasters:
//   the L0 ARM  -- L0's own sampling of the nonlinear v_xc (xcRasterX, booked at L0's x = G_0^2 / p_cap(L0));
//   the RHO ARM -- below the density's rigorous x, the collocated rho aliasing on the FINER levels (rasterX, booked at the finest level's x).
// Both are measured facts.  MnO Gamma uniform cf 4, eps 1e-6 followed L0 ALONE: L0 60^3 -> 72^3 -> 72^3 -> 80^3 (x 9.32, 13.42, 13.42, 16.58)
// read 2.90e-5, 9.98e-6, 9.98e-6 (to 4 digits), 5.39e-6 while the finest level rose 80 -> 96 -> 108 -> 120
// (`~/Code/qchem6-runs/anchors/lane-a-pin10/mno_cf4_pilot.log`).  Na k222s at density x 6 followed rho alone: |dE_xc| 1.2e-4 with the XC
// raster at x 28 AND at x 60, 5e-8 once the density raster went to x 16 (C2 step 3).  One audit cannot tell them apart, so an over XC term
// raises BOTH arms, and a round whose error did not fall REFUTES an arm that moved alone (it stops demanding).  Each arm's demand is the
// power law fitted on ITS OWN axis once that axis has moved (MnO's L0 points fit d ~ x^-2.9 on both intervals; the rigorous e^-x shape,
// which the pilot used alone before, under-asked every round there).  A policy (the pilot's), so it lives in the policy layer.
module;
#include <optional>
#include <string>
export module qchem.RasterPilotXC;

export namespace qchem
{

//! The uniform XC term's two-arm demand across the raster pilot's rounds.  Feed it each round's audit with \c Round, in order.
class UniformXCDemand
{
public:
    //! The error "fell" between two rounds when it dropped below this fraction of the round before; otherwise the arm that moved alone is
    //! refuted.  A judgement threshold, not a calibration: MnO's unmoved L0 re-measured its error to 4 digits (9.981e-6 both rounds) and
    //! every moved round fell by 1.8x or more, so anything between 1.0 and ~0.6 separates them.
    static constexpr double kImprovedRatio=0.9;
    //! At most DOUBLE an arm's x per round: a slope from two nearby, noisy rounds extrapolates wildly (Si: k 1.4 from x 21 -> 30 would ask
    //! x ~ 3000), and the next round measures again anyway.
    static constexpr double kMaxStep=2.0;

    struct Verdict
    {
        double needRho=0.0;      //!< the rho arm's demanded density x (0 = no demand: meets, refuted, or not below the rigorous x)
        double needL0 =0.0;      //!< the L0 arm's demanded L0 x (0 = no demand)
        bool   over   =false;    //!< the delivered error is over the target
        bool   unanswered=false; //!< over, with neither arm left to raise
        bool   rhoArm =false;    //!< the rho arm is in play this round (below the density's rigorous x, not refuted)
        std::string refuted;     //!< the arm refuted THIS round ("density raster" / "L0 raster"), empty when none
    };

    //! One round: \a xRho the finest level's delivered x, \a xL0 L0's, \a delivered the XC term's error (|dE| or max|dH|), \a belowRho
    //! whether \a xRho sits below the density's rigorous x.
    Verdict Round(double xRho, double xL0, double delivered, double target, bool belowRho);

    bool RhoRefuted() const {return itsRhoRefuted;}
    bool L0Refuted () const {return itsL0Refuted;}

private:
    struct Point { double xRho, xL0, d; };
    std::optional<Point> itsBefore;   //!< the round before: each arm's power-law fit's first point
    bool itsRhoRefuted=false, itsL0Refuted=false;
};

} // namespace qchem
