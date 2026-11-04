// File: Calculation/tests/MoleculeSetup.C  The molecular RUN SETUP (doc/CalcNetworkReview.md §7, step 7): the decisions made
// once -- spin promotion, the valence electron count from the pseudopotential set, the parity and accelerator refusals.
#include <gtest/gtest.h>
#include <stdexcept>
#include <nlohmann/json.hpp>
import qchem.MoleculeCalculation;
import qchem.Structure;
import qchem.Types;

using namespace qchem;

namespace
{
Molecule OneAtom(int Z)
{
    Molecule m;
    m.Insert(Site{Z, rvec3_t(0,0,0)});
    return m;
}
} // namespace

TEST(MoleculeSetup, AnOpenShellIsPromotedToPolarizedAndItsParityIsChecked)
{
    const auto s=MoleculeSetup::Resolve(OneAtom(7), {.multiplicity=4});
    EXPECT_EQ(s->opts.spin, SpinGroup::Polarized);
    EXPECT_EQ(s->Ne, 7);
    EXPECT_FALSE(s->pp) << "no pp named: all-electron";
    EXPECT_THROW(MoleculeSetup::Resolve(OneAtom(7), {.multiplicity=3}), std::runtime_error) << "a triplet of 7 electrons";
}

TEST(MoleculeSetup, APseudopotentialSetGivesTheValenceCountAndTheValenceStructure)
{
    const Molecule na=OneAtom(11);
    const auto s=MoleculeSetup::Resolve(na, {.pp="gth-lda-valence"});
    ASSERT_TRUE(s->pp);
    EXPECT_EQ(s->Ne, s->pp->Valence("Na"));
    const auto point=s->PointStructure(na);
    EXPECT_EQ(int(point->GetNumElectrons()), s->Ne);
    EXPECT_EQ(int(s->reference->GetNumElectrons()), 11) << "the reference keeps the TRUE charges; only a point is valence-charged";
}

TEST(MoleculeSetup, AnUnknownAcceleratorIsRefusedAtResolve)
{
    EXPECT_THROW(MoleculeSetup::Resolve(OneAtom(2), {}, {.type="Anderson"}), std::runtime_error);
}

// EVERY KNOB OF THE MOLECULAR CHAIN REACHES THE FACTORY'S DOOR (2026-11-02): the GDM trust radius and the Ladder's hand-off triggers were
// unreachable -- MoleculeCalculation passed only {NProj, EMax, EMin, SVTol}.  ToAcceleratorJson is the one translation both facades use; its
// defaults are the door's own, so a molecule that sets none of them runs what it ran.
TEST(MoleculeSetup, EveryAcceleratorKnobReachesTheFactoryDoor)
{
    const AcceleratorOptions def;
    EXPECT_EQ(ToAcceleratorJson(def, 0.3), (nlohmann::json{{"NProj",4},{"EMax",0.3},{"EMin",1e-7},{"SVTol",5e-9},{"Trust",0.1},
                                                         {"ethresh",1e-8},{"stall",5},{"floor",1e-8},{"switchat",1e-6},{"signal","EnergyChange"}}));
    AcceleratorOptions a; a.nProj=6; a.trust=0.25; a.ladderEThresh=2e-8; a.ladderStall=7; a.ladderFloor=3e-8; a.ladderSwitchAt=1e-6;
    const nlohmann::json j=ToAcceleratorJson(a, 0.5);
    EXPECT_EQ(j.at("NProj"), 6); EXPECT_EQ(j.at("EMax"), 0.5); EXPECT_EQ(j.at("Trust"), 0.25);
    EXPECT_EQ(j.at("ethresh"), 2e-8); EXPECT_EQ(j.at("stall"), 7); EXPECT_EQ(j.at("floor"), 3e-8); EXPECT_EQ(j.at("switchat"), 1e-6);
}
