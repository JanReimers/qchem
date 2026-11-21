// File: Calculation/tests/BasisMenu.C  The orbital-basis menu a structure offers (qchem.BasisMenu; lane G's design note
// laneG-basis-menu-design.md, ruled 2026-11-21): the radial atom's families only for one site, the Gaussian files greyed with the
// reader's reason when they miss an element, the crystal's plane waves.
#include <gtest/gtest.h>
#include <algorithm>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>
import qchem.BasisMenu;
import qchem.Deck;
import qchem.SolidCalculation;   // RunOptions
import qchem.Structure;
import qchem.UnitCell;
import qchem.Types;

using namespace qchem;

namespace
{
const BasisOffer* Find(const std::vector<BasisOffer>& m, const std::string& name)
{
    auto it=std::find_if(m.begin(), m.end(), [&](const BasisOffer& o){ return o.name==name; });
    return it==m.end() ? nullptr : &*it;
}
}

TEST(BasisMenu, OneSiteOffersTheRadialFamiliesAndTheOptInMolecule)
{
    const auto m=OrbitalBasisOffers(Atom(10, 0.0));
    const auto names=RadialAtomBasisNames();
    ASSERT_EQ(m.size(), names.size()+1);
    for (size_t i=0;i<names.size();i++)
    {
        EXPECT_EQ(m[i].route, "atom");
        EXPECT_EQ(m[i].name, names[i]);
        EXPECT_EQ(m[i].key, "atom.type="+names[i]);
        EXPECT_TRUE(m[i].admitted);
        EXPECT_FALSE(m[i].optIn);
        EXPECT_EQ(m[i].note.find("relativistic")!=std::string::npos, names[i].ends_with("_RKB")) << names[i];
    }
    EXPECT_EQ(m.back().route, "molecule");
    EXPECT_TRUE(m.back().optIn) << "the one-atom molecule is offered only when asked (user ruling Q1)";
    EXPECT_TRUE(std::none_of(m.begin(), m.end(), [](const BasisOffer& o){ return o.key.starts_with("basis.data="); }))
        << "a lone atom lists no Gaussian file (the deck refuses a basis block on an atom)";
    Molecule lone; lone.Insert(Site{10, rvec3_t(0,0,0)});
    EXPECT_EQ(OrbitalBasisOffers(lone).size(), m.size()) << "one site is the radial atom, whatever class holds it (deck::KindOf)";
}

TEST(BasisMenu, EachRadialNameIsTheDecksSpelling)
{
    const auto names=RadialAtomBasisNames();
    for (size_t i=0;i<names.size();i++)   // the deck reads each name to the family of the same place in BasisSet::Radial::Type
    {
        RunOptions o;
        ASSERT_NO_THROW(deck::FromJson(nlohmann::json{{"schema",2},{"atom",{{"type",names[i]}}}}, o)) << names[i];
        EXPECT_EQ(o.atom.type, BasisSet::Radial::Type(i)) << names[i];
    }
    RunOptions o;
    EXPECT_ANY_THROW(deck::FromJson(nlohmann::json{{"schema",2},{"atom",{{"type","NoSuchFamily"}}}}, o)) << "the deck refuses an unknown name";
}

TEST(BasisMenu, AMoleculeGreysTheFilesThatMissAnElement)
{
    Molecule water;
    water.Insert(Site{8, rvec3_t(0,0,0)});
    water.Insert(Site{1, rvec3_t(0, 1.431,1.107)});
    water.Insert(Site{1, rvec3_t(0,-1.431,1.107)});
    const auto m=OrbitalBasisOffers(water);
    const BasisOffer* dzvp=Find(m,"dzvp");
    ASSERT_NE(dzvp, nullptr);
    EXPECT_TRUE(dzvp->admitted);
    EXPECT_EQ(dzvp->key, "basis.data=dzvp");
    EXPECT_NE(dzvp->blockGroup.find("C2v"), std::string::npos) << dzvp->blockGroup;
    const BasisOffer* vb=Find(m,"valence_lowq_vb");   // F Na Al Mn O: no H block
    ASSERT_NE(vb, nullptr);
    EXPECT_FALSE(vb->admitted);
    EXPECT_EQ(vb->reason, "valence_lowq_vb has no block for H");
    EXPECT_TRUE(std::all_of(m.begin(), m.end(), [](const BasisOffer& o){ return o.route=="molecule" && !o.optIn; }));
    EXPECT_EQ(MoleculeBasisOffers(water).size(), m.size());
}

TEST(BasisMenu, TheOptInRouteLeadsToTheMoleculesMenu)
{
    const auto m=MoleculeBasisOffers(Atom(10, 0.0));
    const BasisOffer* dzvp=Find(m,"dzvp");
    ASSERT_NE(dzvp, nullptr);
    EXPECT_EQ(dzvp->route, "molecule");
    EXPECT_TRUE(dzvp->admitted) << "dzvp holds Ne";
}

TEST(BasisMenu, ACrystalAddsThePlaneWaves)
{
    FCCUnitCell si(10.26);
    si.AddAtom(14, {0,0,0});
    si.AddAtom(14, {0.25,0.25,0.25});
    const auto m=OrbitalBasisOffers(si);
    ASSERT_FALSE(m.empty());
    const BasisOffer* sr=Find(m,"sipp_sr");
    ASSERT_NE(sr, nullptr);
    EXPECT_TRUE(sr->admitted);
    EXPECT_EQ(m.back().key, "basis.family=planeWave");
    EXPECT_TRUE(std::all_of(m.begin(), m.end(), [](const BasisOffer& o){ return o.route=="crystal" && o.blockGroup=="T (by k)"; }));
}
