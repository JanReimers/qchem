// File: IntegrationTests/GPW/TiO2.C  rutile TiO2 -- the CONVERGENCE PANEL's little-group REFERENCE case (doc/SCFRobustnessPlan.md §4).
//
// littleGroupK measures the ground state's group on the reference orbitals to 1e-6, so its SCF must reach drho/N 1e-9
// (`a2e664f1`, deck::Vet refuses less).  The case is the SCF of `decks/campaigns/TiO2/tio2_k222_q222_littleGroupGate.json` with its
// response post-actions and its state file stripped: can the recipe deliver a 1e-9 reference within budget?  Full tier.
#include "gtest/gtest.h"
#include <fstream>
#include <string>
#include <nlohmann/json.hpp>

import qchem.Tests.Panel;   // panel::Run, panel::BasisDataDir

TEST(GPW_TiO2, k222_Pol_Smear_Kerker_Pulay_Converges_Panel_Long)
{
    std::ifstream in(qchem::panel::DeckDir().parent_path()/"campaigns"/"TiO2"/"tio2_k222_q222_littleGroupGate.json");
    ASSERT_TRUE(in) << "the TiO2 little-group gate deck";
    nlohmann::json deck=nlohmann::json::parse(in);
    deck.erase("postSCF"); deck.erase("state");
    deck["basis"]["file"]=(qchem::panel::BasisDataDir()/"valence_lowq_tm.bsd").string();   // the deck names /home/janr/Code/qchem6/...
    const auto o=qchem::panel::Run(deck, 49, "TiO2 k222 little-group reference to drho/N 1e-9 (regime 1, kT 5 mHa; 39 measured, 2026-10-29 panel baseline)");
    ASSERT_TRUE(o.run.converged) << o.run.summary;
}
