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
