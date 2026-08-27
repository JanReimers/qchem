// File: Calculation/tests/RunDiagnostics.C
//
// THE T2/T3/T4 OUTCOME DETECTORS, AT THE UNIT LEVEL (user, 2026-08-27).
//
// Until now these rules had ZERO unit coverage: grep src/*/tests/ for RunDiagnostics, OrderCollapsed or
// ChargeSloshed and nothing came back.  They were exercised only by coaxing a full SCF into the outcome
// under test -- which is how `GPW_SCF.ImposedOrderLostIsAPostconditionFailure_Na2Box` came to need
// alpha=0.5 and a 100-iteration cap as LOAD-BEARING parameters: not because the physics wants them, but
// because `NotConverged` is checked before `OrderLost`, so the run has to CONVERGE before the branch under
// test can even be reached.
//
// ⇒ THAT MAKES A LOGIC TEST HOSTAGE TO A CONVERGENCE MARGIN -- textbook flakiness (user: "all SCF runs are
// strictly speaking flaky, but this one is much more flaky than average+3sigma").  Measured 2026-08-27: the
// fixture converges in 66 of its 100-iteration cap on the walk and in 290 on the contraction kernel, so a
// change four layers below it silently disarms the gate.  Its own comment already records that this had
// happened once before -- "the postcondition below is unexercised AGAIN".
//
// The rules themselves are PURE FUNCTIONS of two trajectories and two flags, so they belong here, tested
// on synthetic input in microseconds and deterministically.  The integration test's remaining job is
// narrow and honest: prove the WIRING reaches the detector on a real run.
#include "gtest/gtest.h"
#include <vector>
#include <string>

import qchem.SolidCalculation;   // RunDiagnostics

using qchem::RunDiagnostics;

//! The friend declared in src/forward.H: it is what lets a test build a diagnostics object carrying a
//! CHOSEN trajectory instead of whatever an SCF happened to produce.
class RunDiagnosticsTests : public ::testing::Test
{
protected:
    //! Build a RunDiagnostics with an explicit history.
    //! \param order per-iteration integrated site moment, \param seed the RAW seed's moment,
    //! \param eee per-iteration Hartree term, \param converged the run's own verdict.
    static RunDiagnostics Make(std::vector<double> order, double seed, bool hasBasins,
                               std::vector<double> eee, bool converged)
    {
        RunDiagnostics d;
        d.itsOrder     = std::move(order);
        d.itsSeedOrder = seed;
        d.itsHasBasins = hasBasins;
        d.itsEee       = std::move(eee);
        d.itsConverged = converged;
        return d;
    }
};

//============================================================================================
// T2/T4 -- the ORDER channel
//============================================================================================

TEST_F(RunDiagnosticsTests, OrderThatRisesThenDiesIsACollapse)
{
    // The MnO shape the rule was written from: small at first, GROWS as the exchange splitting builds,
    // then decays to nothing.  Judged against iteration 1 this "survived"; against the peak it is a
    // collapse, which is precisely why OrderPeak is a high-water mark.
    const RunDiagnostics d=Make({0.0046,0.02,0.06,0.106,0.04,0.001,7e-5,7e-5}, 0.0, true, {}, true);
    EXPECT_TRUE (d.HasOrder());
    EXPECT_NEAR (d.OrderPeak(), 0.106, 1e-12);
    EXPECT_TRUE (d.OrderCollapsed());
    EXPECT_GT   (d.OrderDiedAt(), 0u) << "a collapse must name the step it died at";
}

TEST_F(RunDiagnosticsTests, TheRawSeedCountsAsAPeakEvenWhenTheFirstFillEatsIt)
{
    // The Na2 defect this rule exists for: the density reachable at construction is ALREADY one aufbau
    // fill downstream, so a seed staggered at +-0.47 e reads ~0.001 e by iteration 1.  Judged on the
    // iterates alone there is no baseline and the collapse is invisible.
    const RunDiagnostics d=Make({0.001,0.0009,0.0008}, 0.4686, true, {}, true);
    EXPECT_NEAR(d.OrderPeak(), 0.4686, 1e-12) << "the RAW seed must be able to be the peak";
    EXPECT_TRUE(d.OrderCollapsed());
}

TEST_F(RunDiagnosticsTests, OrderThatSurvivesIsNotACollapse)
{
    // The Mn2 mirror-image gate: a real magnet that keeps its order.  The detector must DISCRIMINATE,
    // not always fire -- a rule that convicts every run explains none of them.
    const RunDiagnostics d=Make({4.78,4.6,4.4,4.3,4.222}, 4.781, true, {}, true);
    EXPECT_TRUE (d.HasOrder());
    EXPECT_FALSE(d.OrderCollapsed());
    EXPECT_EQ   (d.OrderDiedAt(), 0u) << "a run that never died must not name a step";
}

TEST_F(RunDiagnosticsTests, OrderThatDipsAndRecoversIsNotACollapse)
{
    // "Fell below 1% of the peak" is not enough on its own -- the rule is fell AND STAYED.  An SCF that
    // passes through a near-zero moment on its way to a magnetic answer must not be convicted.
    const RunDiagnostics d=Make({0.5,0.3,0.001,0.2,0.45,0.5}, 0.5, true, {}, true);
    EXPECT_FALSE(d.OrderCollapsed()) << "a transient dip is not a collapse";
    EXPECT_EQ   (d.OrderDiedAt(), 0u);
}

TEST_F(RunDiagnosticsTests, NoSiteBasinsMeansTheCheckSkipsRatherThanGuesses)
{
    // A uniform XC mesh owns no atom-centred basins, so SiteMoments returns empty and the INTEGRATED
    // moment is not measurable at all.  The detector must decline, never read the absence as m=0 --
    // that would convict every uniform-mesh run of losing order it never had.
    const RunDiagnostics d=Make({}, 0.0, false, {}, true);
    EXPECT_FALSE(d.HasOrder());
    EXPECT_FALSE(d.OrderCollapsed()) << "no basins must SKIP, not fire";
}

TEST_F(RunDiagnosticsTests, AnUnpolarisedRunHasNoOrderToLose)
{
    const RunDiagnostics d=Make({0.0,0.0,0.0}, 0.0, true, {}, true);
    EXPECT_FALSE(d.HasOrder());
    EXPECT_FALSE(d.OrderCollapsed());
}

//============================================================================================
// T3 -- the CHARGE channel
//============================================================================================

TEST_F(RunDiagnosticsTests, ChargeSloshIsJudgedOnTheEndNotThePeak)
{
    // The calibration that mattered (doc/OpenWork.md T3): a healthy run can RESTRUCTURE on its way to
    // the answer and raise Eee a long way before coming back.  The MnO free run peaked at 39.98 and
    // returned to 15.64 -- a peak-based rule reads 3.15x and convicts it; the END-based rule reads 1.23x
    // and correctly stays silent.
    const RunDiagnostics d=Make({}, 0.0, false, {13.48,20.0,39.98,25.0,15.64}, false);
    EXPECT_LT  (d.SloshRatio(), RunDiagnostics::kSloshFactor);
    EXPECT_FALSE(d.ChargeSloshed()) << "a run that came back down is not sloshing";
}

TEST_F(RunDiagnosticsTests, ChargeThatRunsAwayAndStaysThereIsSlosh)
{
    const RunDiagnostics d=Make({}, 0.0, false, {13.48,20.0,29.0,35.1}, false);
    EXPECT_GT  (d.SloshRatio(), RunDiagnostics::kSloshFactor);
    EXPECT_TRUE(d.ChargeSloshed());
}

TEST_F(RunDiagnosticsTests, AConvergedRunNeverSloshesByConstruction)
{
    // The containment that makes a mistake in this detector cost a LABEL and never a good answer: a
    // converged density is stationary, so the detector is confined to runs that already failed.
    const RunDiagnostics d=Make({}, 0.0, false, {13.48,20.0,29.0,35.1}, /*converged*/true);
    EXPECT_FALSE(d.ChargeSloshed()) << "ChargeSloshed must be false on a converged run by construction";
}
