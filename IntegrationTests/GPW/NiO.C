// File: IntegrationTests/GPW/NiO.C  NiO AFM-II -- the CONVERGENCE PANEL's restored-stage recipe case (doc/SCFRobustnessPlan.md §4, Q-R).
//
// Q-R ruling (b) (user 2026-11-13): a NAMED MagneticInsulator (or Metal) recipe RESTORES on its own Null + Kerker+Pulay history from
// iteration 1 (pulayStart 0), not on row 4's DIIS + linear D-mixing <= 0.3 (src/Calculation/SCFRecipe.C).  The case is a WARM restart:
// a vendored pilot-sized state restored onto other grids, so the density MOVES.  An EXACT resume of the same state converges in 1
// iteration under either recipe (measured 2026-11-13), so it cannot tell them apart and is not the gate.
#include "gtest/gtest.h"
#include <cstdio>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>

import qchem.Tests.Panel;   // panel::LoadDeck (decks/panel/, a relative restartFrom resolved against decks/states/)
import qchem.Deck;          // deck::Build / deck::ScheduleToRun (the observer sees every iterate)
import qchem.SCFCap;        // IterationBudget
import qchem.SCFIterator;   // IterationRow (the first-iterate probe)
import qchem.SolidCalculation;
import qchem.Reporting;
import qchem.Types;

// Vendored state decks/states/NiO_AFM2_k222_MagneticInsulator.h5: saved by ~/Code/qchem6-runs/NiO/laneA_qr_cold.json (MagneticInsulator
// cold, 18 iterations on the ball rasters, then the raster pilot's 10-iteration re-convergence onto cutoffFactor 4 / V_long x 4.97,
// E -106.3510571976), 2026-11-13.  The deck restores it onto the cold run's FIRST grids (cutoffFactor 2, rasters x 0), whose converged
// energy is that same cold run's -106.3514267289 (Ni moment +-1.376): the independent route for the fixed point below.
// MEASURED 2026-11-13 (6 threads): 11 iterations, first iterate -106.351423122 (3.6 uHa above), E -106.3514267282.  The SAME restart as
// recipe Insulator (row 4, DIIS + linear 0.3) stops at iteration 1 on a drho of 6.5e-13 while [F,D] is 6.9e-4, 147 nHa high (OpenWork §3
// "A Δρ-only stop certifies an unconverged iterate") -- the case for ruling (b).
// The state's OWNER is decks/panel/NiO_k222_MagneticInsulator_Cold.json (~1.5 h, run by hand, not a gtest).  REFRESH IT ONLY WHEN THIS GATE'S
// CLAIM REQUIRES IT (controller ruling 2026-11-13), never as routine: every refresh adds the full 6.4 MB to git history.
TEST(GPW_NiO, k222_Pol_Smear_Kerker_Pulay_Null_Saved_RestoresOnItsOwnHistory_Panel_Long)
{
    qchem::deck::RunSpec spec;
    qchem::deck::FromJson(qchem::panel::LoadDeck("NiO_k222_MagneticInsulator_WarmRestore.json"), spec);
    const qchem::report::Run run("NiO_k222_warmRestore");
    std::vector<double> E;
    auto calc=qchem::deck::Build(spec, [&E](const qchem::SCFIterator::IterationRow& p){ E.push_back(p.eb.GetTotalEnergy()); });
    ASSERT_TRUE(calc);
    std::size_t used=0;
    bool converged=false;
    {
        qchem::IterationBudget b(13, "NiO k222 MagneticInsulator warm restore onto the cold grids (11 measured, 2026-11-13)");
        converged=bool(calc->Converge(qchem::deck::ScheduleToRun(spec)));
        used=b.Used();
    }
    ASSERT_TRUE(converged);
    ASSERT_FALSE(E.empty());
    std::printf("[Q-R] iterations %zu  first iterate %.12f  E %.12f\n", used, E.front(), calc->Result()->Energy());
    // The FIRST iterate (CLAUDE.md: pin a restart's first iterate): the restored D's energy on these grids, before any mixing.
    EXPECT_NEAR(E.front(), -106.351423122, 1e-8) << "the first iterate of the warm restore";
    EXPECT_GE(used, 2u) << "a moved density cannot be certified at iteration 1 (row 4 stopped there, 147 nHa high)";
    EXPECT_NEAR(calc->Result()->Energy(), -106.3514267289, 1e-8) << "the cold run's fixed point on these grids";
    const auto mu=calc->LastIterateTerms().charge.siteMoments;
    ASSERT_EQ(mu.size(), 4u) << "one Becke basin per site (Ni, Ni, O, O)";
    EXPECT_NEAR(mu[0],  1.376, 0.01) << "AFM-II survives the restore";
    EXPECT_NEAR(mu[1], -1.376, 0.01);
}
