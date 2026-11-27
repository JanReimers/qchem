// File: IntegrationTests/GPW/MnO2.C  beta-MnO2 A1AFM -- the CONVERGENCE PANEL's regime-2'' cases (doc/SCFRobustnessPlan.md §4).
//
// U on the O p_z (plane-normal, slot Gamma2) direction softens A1AFM toward the ferrimagnetic state 11.7 mHa below it at U = 0
// (OpenWork §3, QE agrees), so the small-U leg is where the density mixers fail: Kerker+Pulay overshoots the U = 1 eV basin and
// plain Kerker limit-cycles at U = 2 eV (OpenWork §3 "Kerker+Pulay overshoots β-MnO₂'s U = 1 eV basin").  Each case is a
// vendored deck in decks/panel/ run COLD from IonicSAD, with its cap at the target budget: a red case fails fast.
// Full tier (_Panel_Long): ~20 s/iteration at 8 threads.
#include "gtest/gtest.h"
#include <string>
#include <nlohmann/json.hpp>

import qchem.Tests.Panel;   // panel::LoadDeck / panel::Run (decks/panel/, under an IterationBudget)

namespace {
// The control-basis U = 1 eV state, plain Kerker 0.1 warm from U = 0.5 (`~/Code/qchem6-runs/MnO2/uscan/control_U1a.log`, r018,
// 2026-10-29): the fixed point both U = 1 cases must land on.
constexpr double kE_U1=-93.36627498;
}

// RED (QCHEM_PANEL_RED): Kerker+Pulay (history 8, alpha 0.45) from the seed reaches -93.3663 at drho 4.5e-4 and then runs away
// (E UP 0.1 Ha, repeatedly; `uscan/prefix_restartbug/control_U1a.pulay8.log`).  Target: U = 0 / 4.62 eV converge in 16-24.
TEST(GPW_MnO2, Γ_U_Shub_Pol_Smear_Kerker_Pulay_Converges_Panel_Long)
{
    const auto o=qchem::panel::Run(qchem::panel::LoadDeck("MnO2_U1_KerkerPulay.json"), 40,
                                   "beta-MnO2 U(O p_z) = 1 eV, Kerker+Pulay, cold (regime 2''; OpenWork §3 overshoot row)");
    ASSERT_TRUE(o.run.converged) << o.run.summary;
    ASSERT_TRUE(o.run.energy);
    EXPECT_NEAR(*o.run.energy, kE_U1, 1e-4) << "the U = 1 eV A1AFM basin";
}

// Plain Kerker alpha 0.1 held the U = 1 eV basin WARM; this is the COLD start (143 measured, 2026-10-29 panel baseline; the +U pseudo-atom solves are setup, not counted).
TEST(GPW_MnO2, Γ_U_Shub_Pol_Smear_Kerker_Converges_Panel_Long)
{
    const auto o=qchem::panel::Run(qchem::panel::LoadDeck("MnO2_U1_Kerker.json"), 179,
                                   "beta-MnO2 U(O p_z) = 1 eV, plain Kerker 0.1, cold (regime 2''; 143 measured, 2026-10-29 panel baseline)");
    ASSERT_TRUE(o.run.converged) << o.run.summary;
    ASSERT_TRUE(o.run.energy);
    EXPECT_NEAR(*o.run.energy, kE_U1, 1e-4) << "the U = 1 eV A1AFM basin";
}

// RED: plain Kerker 0.1 at U = 2 eV limit-cycles (E -93.336 ... -93.352, cell moment 3.35-3.44, 300 iterations, warm;
// `uscan/control_U2a.log`); the p5f3_Od basis converged there in 78.  No fixed point is banked on this basis: convergence only.
TEST(GPW_MnO2, Γ_U_Shub_Pol_Smear_Kerker_ConvergesAtU2eV_Panel_Long)
{
    const auto o=qchem::panel::Run(qchem::panel::LoadDeck("MnO2_U2_Kerker.json"), 80,
                                   "beta-MnO2 U(O p_z) = 2 eV, plain Kerker 0.1, cold (regime 2''; OpenWork §3 overshoot row)");
    ASSERT_TRUE(o.run.converged) << o.run.summary;
}

// PANEL case 19 (regime 4, EXACT RESUME): the control-basis U = 0 state r015, vendored (decks/states/), resumed on its own recipe.  Iteration 1 must
// open AT the saved state and the resume must not wander: before D5 (2026-10-29, the carried D was not the working density) iteration 1 opened
// 4 mHa off and the resume took 23 iterations; fixed: 4.
TEST(GPW_MnO2, Γ_U_Shub_Pol_Smear_Kerker_Pulay_Saved_eqExactResume_Panel_Long)
{
    const auto o=qchem::panel::Run(qchem::panel::LoadDeck("MnO2_U0_ExactResume.json"), 5,
                                   "beta-MnO2 U = 0 exact resume of the vendored r015 state (regime 4; 4 measured, 2026-10-29)");
    ASSERT_TRUE(o.run.converged) << o.run.summary;
    ASSERT_TRUE(o.run.energy);
    EXPECT_NEAR(*o.run.energy, -93.3837082, 1e-6) << "an exact resume lands on the state it resumed";
}

// The DENSE-d basis (lane H's p5f3_dD, Mn d midpoints 0.56 ... 24.66, vendored): S's smallest eigenvalue is 6.75e-7, and the orthonormal frame
// amplifies a raster alias by up to 1/lambda_min.  On main 0fab43d6 both decks below DIVERGED from iteration 1 (+4.2e4 / +845 Ha) and ran 300
// iterations to +4.7e3 Ha (OpenWork_History5 "Dense Mn d midpoints").  Measured 2026-11-01 (lane A), first iterate vs raster x: 8 -> +41888,
// 11 -> -80.75302, 13 / 16 -> -80.75311.  A periodic first iterate above +1 Ha per electron is refused, not iterated (here +1102 / +22).
// RE-MEASURED 2026-11-24 (lane A, the absolute pair routing ecut >= 2 p ln(1/eps)): x 8 no longer aliases (first iterate -80.75, converges in
// 25), so the refusal case states x 4, where the first iterate is still +2080.3 Ha (54.74 per electron; x 2 the same).

// A STATED raster is the deck's choice: never overridden, refused at iteration 1 with the cause named (~3 min, mostly setup).
TEST(GPW_MnO2, Γ_U_Shub_Pol_Smear_Kerker_Pulay_RefusesAnAliasedFirstIterate_Long)
{
    const auto o=qchem::panel::Run(qchem::panel::LoadDeck("MnO2_dD_StatedX4.json"), 1,
                                   "beta-MnO2 dense-d at a stated x 4: the first iterate is refused, not iterated");
    EXPECT_FALSE(o.run.converged);
    EXPECT_FALSE(o.run.energy) << "a refused run has no converged energy";
    EXPECT_EQ(o.iterations, 1u) << "refused AT iteration 1";
    EXPECT_NE(o.run.summary.find("per electron"), std::string::npos) << o.run.summary;
}

// At the DEFAULT accuracy the first stage starts on the ball's own rasters (x 4: Etot +845 Ha), is refused, and the pilot restarts ONCE on
// rasters whose alias does not exceed lambda_min(S) (x 14.2, 90^3).  With the density ball judged (2026-11-03) the pilot doubles the ball once
// (cf 2 -> 4, Eee 3.3e-4 OVER 6e-5) and re-converges: 47 measured = 25 + 22 re-convergence (was 25 with 0); E moves +334 uHa (ball).
// RE-JUDGED 2026-11-25 (doc/Pins.md pin 10; lane B, D-BECKE-ENERGY-AUDIT fallout): the energy run's pilot now climbs the Becke XC mesh
// (40,L29) -> (80,L35) and E moves +858.5 uHa, -93.4416050775 -> -93.4407465622.  Two independent routes, same deck, same rasters (x 14.17):
//   - Becke STATED at (40,L29): -93.44160394, 1.1 uHa from the old anchor, and the pilot names its XC term OVER the 6e-5 target (NOT SETTLED).
//     So the old anchor was the entry mesh's quadrature error, not the state.
//   - the UNIFORM XC route (xc.mesh Uniform): -93.44073782, SETTLED, 8.7 uHa from the climbed Becke value.
// The routes agree well inside the target; the climbed value is re-banked.  Decks and logs: ~/Code/qchem6-runs/MnO2/laneB_pin10/.
TEST(GPW_MnO2, Γ_U_Shub_Pol_Smear_Kerker_Pulay_ConvergesThroughARefusedFirstIterate_Panel_Long)
{
    const auto o=qchem::panel::Run(qchem::panel::LoadDeck("MnO2_dD_Default.json"), 59,
                                   "beta-MnO2 dense-d, default accuracy: refused on the ball rasters, restarted on conditioned ones, one ball doubling (47 measured, 2026-11-03)");
    ASSERT_TRUE(o.run.converged) << o.run.summary;
    ASSERT_TRUE(o.run.energy);
    EXPECT_NEAR(*o.run.energy, -93.4407465622, 1e-5) << "the dense-d A1AFM U = 0 state (re-judged 2026-11-25 on the climbed Becke mesh, uniform route 8.7 uHa away; was -93.4416050775 on (40,L29))";
}
