// File: Calculation/RasterPilotXC.C  The raster pilot's demand rule for the UNIFORM XC term (lane A, 2026-11-25/26; OpenWork "the raster
// pilot never weighs a bigger ball against a finer raster, and its XC demand underestimates every jump").
//
// v_xc[rho] on the uniform route is sampled on ladder level L0, so the XC term's QUADRATURE error is L0's: xcRasterX, booked at L0's
// x = G_0^2 / p_cap(L0).  That is the pilot's ONE XC knob (user ruling 2026-11-26).  The density's INTERPOLATION onto L0's points is the
// SAMPLER's to certify (pin 28: the routing bound, SolidCalculation::AuditXCSampler), never a pilot knob.  Measured 2026-11-26 at a fixed
// L0, the density raster moves the converged E by 34 nHa on MnO, 5 nHa on Si Gamma, 1.2 uHa on Na k222s: the old "XC via rho" arm was
// L0 in disguise (L0 is built at max(rasterX, xcRasterX)) or the sampler's own miss (NaF / Al below density x 8).
//
// The arm's demand is the POWER LAW fitted on L0's own axis once L0 has moved (MnO Gamma uniform cf 4, eps 1e-6: L0 60^3 -> 72^3 -> 80^3,
// x 9.32, 13.42, 16.58, read 2.90e-5, 9.98e-6, 5.39e-6 -- d ~ x^-2.9 on both intervals; `~/Code/qchem6-runs/anchors/lane-a-pin10/
// mno_cf4_pilot.log`), the rigorous e^-x shape before that (it under-asks every round on its own).  A round in which L0 MOVED and the
// error did not fall REFUTES the arm: raising L0 does not answer this error, and the term is STALLED (lane B's selector A climbs to a Becke
// mesh on it).  A round in which L0 did not move (FFT rounding) refutes nothing.  A policy (the pilot's), so it lives in the policy layer.
module;
#include <optional>
export module qchem.RasterPilotXC;

export namespace qchem
{

//! The uniform XC term's L0 demand across the raster pilot's rounds.  Feed it each round's audit with \c Round, in order.
class UniformXCDemand
{
public:
    //! The error "fell" between two rounds when it dropped below this fraction of the round before; otherwise a moved L0 is refuted.  A
    //! judgement threshold, not a calibration: MnO's unmoved L0 re-measured its error to 4 digits (9.981e-6 both rounds) and every moved
    //! round fell by 1.8x or more, so anything between 1.0 and ~0.6 separates them.
    static constexpr double kImprovedRatio=0.9;
    //! At most DOUBLE L0's x per round: a slope from two nearby, noisy rounds extrapolates wildly (Si: k 1.4 from x 21 -> 30 would ask
    //! x ~ 3000), and the next round measures again anyway.
    static constexpr double kMaxStep=2.0;

    struct Verdict
    {
        double needL0 =0.0;      //!< the demanded L0 x (0 = no demand: meets, or refuted)
        bool   over   =false;    //!< the delivered error is over the target
        bool   refuted=false;    //!< L0 was refuted THIS round (it moved, and the error did not fall)
    };

    //! One round: \a xL0 L0's delivered x, \a delivered the XC term's error (|dE| or max|dH|) against \a target.
    Verdict Round(double xL0, double delivered, double target);

    bool L0Refuted() const {return itsRefuted;}
    //! \brief STALLED: the last round audited was over its target and the L0 arm is refuted -- raising L0 did not answer this XC error, so
    //! the uniform route has no knob left.  The ONE name a route switch reads (lane B's selector A climbs to a Becke mesh on it), never a
    //! raw "the error did not fall", which a round whose L0 did not move also shows (MnO, NaF).  A stable public answer of the demand.
    bool Stalled() const {return itsStalled;}

private:
    struct Point { double xL0, d; };
    std::optional<Point> itsBefore;   //!< the round before: the power-law fit's first point
    bool itsRefuted=false;
    bool itsStalled=false;
};

} // namespace qchem
