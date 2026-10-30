// File: Calculation/tests/MoleculeSetup.C  The molecular RUN SETUP (doc/CalcNetworkReview.md §7, step 7): the decisions made
// once -- spin promotion, the valence electron count from the pseudopotential set, the parity and accelerator refusals.
#include <gtest/gtest.h>
#include <stdexcept>
import qchem.Calculation;
import qchem.Structure;
import qchem.Types;

using namespace qchem;

namespace
{
Molecule OneAtom(int Z)
{
    Molecule m;
    m.Insert(new Atom(Z, rvec3_t(0,0,0)));
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
