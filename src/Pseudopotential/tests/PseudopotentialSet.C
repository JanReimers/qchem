// File: src/Pseudopotential/tests/PseudopotentialSet.C  The named pseudopotential sets (pseudopotentials.json): the contract of the data file.
#include "gtest/gtest.h"
#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

import qchem.Pseudopotential.Set;
import qchem.Pseudopotential.GTH_Potentials;
import qchem.PeriodicTable;

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
    try { s.Valence("Fr"); FAIL() << "expected a throw"; }   // the GTH database has no Fr: neither an entry nor the rule covers it
    catch (const std::runtime_error& e) { EXPECT_NE(std::string(e.what()).find("gth-lda-valence"), std::string::npos); }
    try { PseudopotentialSet::Get("nope"); FAIL() << "expected a throw"; }
    catch (const std::runtime_error& e) { EXPECT_NE(std::string(e.what()).find("known:"), std::string::npos); }
}

TEST(PseudopotentialSet, WithReplacesOrAddsAndKeepsNameAndFunctional)
{
    const PseudopotentialSet base=PseudopotentialSet::Get("gth-lda-valence");
    const PseudopotentialSet s=base.With({{"Mn",15},{"Fr",9}});
    EXPECT_EQ(s.Valence("Mn"), 15);  EXPECT_EQ(s.Valence("Fr"), 9);   // replaced / added (an explicit entry wins even where the database has none)
    EXPECT_EQ(s.Valence("O"), base.Valence("O"));                     // the rest untouched
    EXPECT_EQ(s.Name(), base.Name());  EXPECT_EQ(s.Functional(), base.Functional());
    EXPECT_EQ(base.Valence("Mn"), 7) << "With returns a new set; the original is not modified";
}

// THE RULE (OpenWork §3 "A deck cannot name Te with a PP", 2026-11-11): a set resolves EVERY element the GTH database covers at its functional,
// to the smallest q, and each resolved potential carries the functional it was fitted with.  Over the whole database, never a count.
TEST(PseudopotentialSet, EveryElementTheDatabaseCoversResolvesToItsSmallestValenceAndCarriesItsFittedFunctional)
{
    const PseudopotentialSet s=PseudopotentialSet::Get("gth-lda-valence");
    size_t seen=0;
    for (size_t Z=1; Z<=thePeriodicTable().GetNumElements(); Z++)
    {
        const std::string el=thePeriodicTable().GetSymbol(Z);
        const auto qs=GetGTHValences(el, s.Functional());
        EXPECT_EQ(s.Covers(el), !qs.empty()) << el;
        if (qs.empty()) continue;
        ++seen;
        EXPECT_EQ(s.Valence(el), qs.front()) << el;
        EXPECT_EQ(GetGTH(el, s.Functional(), s.Valence(el)).zion, s.Valence(el)) << el;
        EXPECT_EQ(s.FittedFunctional(el), "PADE") << el << ": every GTH LDA entry is CP2K's GTH-PADE";
    }
    EXPECT_GT(seen, 0u);
    EXPECT_EQ(s.Valence("Te"), 6) << "the row's case: Te resolves on defaults";
    for (const std::string f : {"PBE","BLYP"})   // a non-LDA family: the fitted functional is the lookup key itself
        for (const auto q : GetGTHValences("Si", f)) EXPECT_EQ(GetGTH("Si", f, q).fitted, f) << f << " q" << q;
}

// The valences every banked anchor ran at (the old hand-written list) are unchanged by the rule: bit-identical potentials.
TEST(PseudopotentialSet, TheBankedValencesAreTheRulesValences)
{
    const PseudopotentialSet s=PseudopotentialSet::Get("gth-lda-valence");
    const std::vector<std::pair<std::string,int>> banked={{"Al",3},{"C",4},{"Cs",1},{"F",7},{"Fe",8},{"H",1},{"I",7},{"Li",1},{"Mn",7},
                                                          {"N",5},{"Na",1},{"Ni",10},{"O",6},{"S",6},{"Si",4},{"Sr",2},{"Ti",4},{"V",5}};
    for (const auto& [el,q] : banked) EXPECT_EQ(s.Valence(el), q) << el;
}
