// File: Calculation/tests/AtomSetup.C  The radial atom's RUN SETUP (V-FACADE-PARITY P1, user ruling R3): the decisions made once -- the
// atom's view of the options, the accelerator by name (an unknown one refused, the `directmin` alias kept), the seed and its key.
#include <gtest/gtest.h>
#include <stdexcept>
#include <nlohmann/json.hpp>
import qchem.AtomCalculation;
import qchem.Structure;
import qchem.Types;
import qchem.ChargeDensity.Seed;   // SeedStrategy

using namespace qchem;
using AType=SCFAccelerators::Type;

TEST(AtomSetup, TheStructureGivesZAndTheElectronCount)
{
    const auto s=AtomSetup::Resolve(Atom(26, 1.0));
    EXPECT_EQ(s->Z, 26);
    EXPECT_EQ(s->Ne, 25) << "Fe+";
    EXPECT_DOUBLE_EQ(s->structure->GetNetCharge(), 1.0);
    EXPECT_EQ(*s->opts.model, Model::HF) << "the atom's default model";
    EXPECT_EQ(s->opts.seed, ChargeDensity::SeedStrategy::CoreGuess) << "an atom seeds from the core guess unless told otherwise";
    EXPECT_EQ(s->accelerator, AType::DIIS);
    EXPECT_FALSE(s->Polarized());
    const auto q=AtomSetup::Resolve(26, 1);
    EXPECT_EQ(q->Ne, s->Ne) << "the (Z, charge) door builds the same atom";
}

TEST(AtomSetup, TheAcceleratorIsResolvedByNameAndAnUnknownOneIsRefused)
{
    RunOptions o;
    o.atom.accelerator={{"type","Ladder"}};
    EXPECT_EQ(AtomSetup::Resolve(8, 0, o)->accelerator, AType::Ladder) << "the json escape hatch names it";
    RunOptions gdm; gdm.accelerator=AType::GDM;
    EXPECT_EQ(AtomSetup::Resolve(8, 0, gdm)->accelerator, AType::GDM) << "or the options do (its one home)";

    o.atom.accelerator={{"type","directmin"}};
    const auto dm=AtomSetup::Resolve(8, 0, o);
    EXPECT_EQ(dm->accelerator, AType::GDM);
    EXPECT_TRUE(dm->directMin) << "the alias is kept, as GDM from the first iteration";

    o.atom.accelerator={{"type","Bogus"}};
    const auto bogus=AtomSetup::Admit(Atom(8), o);
    ASSERT_FALSE(bogus) << "an unknown name is refused (it used to run DIIS silently)";
    EXPECT_EQ(bogus.Error().policy, "AcceleratorType");

    o.atom.accelerator={{"type","Ladder"}};
    o.accelerator=AType::GDM;
    const auto both=AtomSetup::Admit(Atom(8), o);
    ASSERT_FALSE(both) << "the options and the escape hatch naming different accelerators";
    EXPECT_EQ(both.Error().policy, "AcceleratorType");
}

TEST(AtomSetup, ThePseudoAtomKeysTheSeedLibraryAtItsValence)
{
    const auto s=AtomSetup::Resolve(14, 10, {.atom={.pseudopotential=true}});
    EXPECT_EQ(s->Ne, 4);
    EXPECT_EQ(s->PseudoValence(), 4);
    EXPECT_EQ(s->seedKey.functional, "LDA");
    EXPECT_EQ(s->seedKey.Q(14), 4) << "the pseudo-atom at its valence";
    EXPECT_EQ(AtomSetup::Resolve(14)->seedKey.Q(14), -1) << "all-electron: unkeyed";
}

TEST(AtomSetup, ARadialOnlyMeshSpelledTheOldWayReadsAsWhatItIs)
{
    RunOptions o; o.xcMesh=kAtomMeshDefault; o.xcMesh->angular=qcMesh::AngularKind::Lebedev; o.xcMesh->angularDegree=0;
    const auto s=AtomSetup::Resolve(8, 0, o);
    EXPECT_EQ(s->xcMesh.angular, qcMesh::AngularKind::None);
    EXPECT_EQ(s->opts.xcMesh->angular, qcMesh::AngularKind::None);
}

TEST(AtomSetup, OnlyOneAtomIsAdmitted)
{
    Molecule two;
    two.Insert(Site{1, rvec3_t(0,0,0)});
    two.Insert(Site{1, rvec3_t(0,0,1.4)});
    const auto r=AtomSetup::Admit(two);
    ASSERT_FALSE(r);
    EXPECT_EQ(r.Error().policy, "AtomStructure");
    EXPECT_THROW((void)AtomSetup::Resolve(two), std::invalid_argument) << "the code door throws it";
}

// The radial atom's Coulomb fit is its own exponent pool's (the fit-rung ladder), so a Gaussian fit FILE is refused by name, never ignored.
TEST(AtomSetup, ACoulombFitFileIsRefused)
{
    RunOptions o; o.model=Model::LDA; o.coulomb.auxBasis="A1_coul";
    const auto r=AtomSetup::Admit(Atom(8), o);
    ASSERT_FALSE(r);
    EXPECT_EQ(r.Error().policy, "ForeignOption");
}
