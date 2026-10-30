// File: src/Pseudopotential/tests/PseudopotentialSet.C  The named pseudopotential sets (pseudopotentials.json): the contract of the data file.
#include "gtest/gtest.h"
#include <algorithm>
#include <stdexcept>
#include <string>

import qchem.Pseudopotential.Set;
import qchem.Pseudopotential.GTH_Potentials;

using namespace qchem::Pseudopotential;

// EVERY set, EVERY element: the variant it names exists in the GTH database for the set's functional, and its Zion is that q.
// (A set pointing at a variant the database lacks would otherwise fail only when a run reaches that element.)
TEST(PseudopotentialSet, EveryEntryNamesAnExistingGTHVariantWithThatValence)
{
    const auto names=PseudopotentialSet::Names();
    EXPECT_FALSE(names.empty());
    for (const auto& n : names)
    {
        const PseudopotentialSet s=PseudopotentialSet::Get(n);
        EXPECT_EQ(s.Name(), n);
        EXPECT_FALSE(s.Functional().empty()) << n;
        for (const std::string el : {"Al","Cs","F","Fe","I","Li","Mn","Na","Ni","O","S","Si","Sr","Ti","V"})   // the elements the shipped materials use
        {
            if (!s.Covers(el)) continue;
            const int q=s.Valence(el);
            const auto have=GetGTHValences(el, s.Functional());
            EXPECT_TRUE(std::find(have.begin(), have.end(), q)!=have.end()) << n << ' ' << el << " q" << q << " is not in the database";
            EXPECT_EQ(GetGTH(el, s.Functional(), q).zion, q) << n << ' ' << el;
        }
    }
}

TEST(PseudopotentialSet, SpeciesFollowTheRequestedOrderAndAnUncoveredElementThrowsNamingTheSet)
{
    const PseudopotentialSet s=PseudopotentialSet::Get("gth-lda-valence");
    const auto sp=s.Species({"O","Mn"});
    ASSERT_EQ(sp.size(), 2u);
    EXPECT_EQ(sp[0].first, "O");  EXPECT_EQ(sp[0].second, 6);
    EXPECT_EQ(sp[1].first, "Mn"); EXPECT_EQ(sp[1].second, 7);
    try { s.Valence("Xe"); FAIL() << "expected a throw"; }
    catch (const std::runtime_error& e) { EXPECT_NE(std::string(e.what()).find("gth-lda-valence"), std::string::npos); }
    try { PseudopotentialSet::Get("nope"); FAIL() << "expected a throw"; }
    catch (const std::runtime_error& e) { EXPECT_NE(std::string(e.what()).find("known:"), std::string::npos); }
}
