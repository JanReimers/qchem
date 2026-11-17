// File: Calculation/tests/SetupDecisions.C  The structure-neutral setup decisions (doc/OOD-SOLID-Cleanup.md V-FACADE-SETUP): ONE rule
// every setup calls, so a molecule and a solid answer the same (multiplicity, electron count) the same way.
#include <gtest/gtest.h>
#include <optional>
#include <stdexcept>
#include <utility>
import qchem.SetupDecisions;
import qchem.MoleculeCalculation;   // MoleculeSetup, CalcOptions
import qchem.SolidSetup;            // SolidCalcOptions, SpinChannelCounts
import qchem.Structure;
import qchem.Types;
import qchem.SCFAccelerator.Factory;

using namespace qchem;

//! The rule itself, over every small (multiplicity, Ne): minimal spin, the singlet, and a stated 2S, with the refusals where they belong.
TEST(SetupDecisions, TheSpinRuleOverEverySmallCount)
{
    for (int nel=0; nel<=9; nel++)
        for (int m=-1; m<=6; m++)
        {
            const int twoS = m>1 ? m-1 : nel%2;
            const bool possible = !(m==1 && nel%2!=0) && twoS<=nel && (nel-twoS)%2==0;
            if (!possible)
            {
                EXPECT_THROW((void)ResolveSpinChannels(m, nel, std::nullopt, "test"), std::invalid_argument) << "m " << m << ", Ne " << nel;
                continue;
            }
            const SpinChannels c=ResolveSpinChannels(m, nel, std::nullopt, "test");
            EXPECT_EQ(c.nUp+c.nDown, nel);
            EXPECT_EQ(c.nUp-c.nDown, m==1 ? 0 : twoS) << "m " << m << ", Ne " << nel;
        }
    // a shared Fermi level (a metal) lets a singlet stand on an odd count: the counts are only the starting split
    const SpinChannels metal=ResolveSpinChannels(1, 3, true, "test");
    EXPECT_EQ(metal.nUp, 2); EXPECT_EQ(metal.nDown, 1);
}

//! THE PARITY THE ROW WAS ABOUT: the solid's and the molecule's spin bookkeeping are the SAME rule, refusals included.
TEST(SetupDecisions, TheMoleculeAndTheSolidCountSpinTheSameWay)
{
    for (int Z : {3, 7, 8})   // Li, N, O: odd and even counts
        for (int m : {0, 1, 2, 3, 4})
        {
            Molecule atom; atom.Insert(Site{Z, rvec3_t(0,0,0)});
            const auto mol=MoleculeSetup::Admit(atom, {.basis="dzvp", .multiplicity=m});
            bool solidOK=true; std::pair<int,int> solid{0,0};
            try { solid=SpinChannelCounts(SolidCalcOptions{.multiplicity=m}, Z); } catch (const std::invalid_argument&) { solidOK=false; }
            ASSERT_EQ(bool(mol), solidOK) << "Z " << Z << ", m " << m;
            if (!mol) { EXPECT_EQ(mol.Error().policy, "SpinMultiplicity"); continue; }
            EXPECT_EQ(mol.Value()->nUp,   solid.first)  << "Z " << Z << ", m " << m;
            EXPECT_EQ(mol.Value()->nDown, solid.second) << "Z " << Z << ", m " << m;
        }
}

TEST(SetupDecisions, AnAcceleratorIsStatedByOneOfThreeNames)
{
    using T=SCFAccelerators::Type;
    EXPECT_EQ(AcceleratorTypeByName("DIIS",   "test"), T::DIIS);
    EXPECT_EQ(AcceleratorTypeByName("GDM",    "test"), T::GDM);
    EXPECT_EQ(AcceleratorTypeByName("Ladder", "test"), T::Ladder);
    EXPECT_THROW((void)AcceleratorTypeByName("Pulay", "test"), std::invalid_argument) << "never a silent DIIS";
}

TEST(SetupDecisions, TheDIISEngageGate)
{
    EXPECT_EQ(DIISEngageGate(0.5, 26, true),  0.5) << "a stated gate is honoured";
    EXPECT_EQ(DIISEngageGate(0.0, 26, true),  100.0) << "a density Fock engages from the start";
    EXPECT_DOUBLE_EQ(DIISEngageGate(0.0, 10, false), 100*0.08/32) << "the Z-scaled gate";
}
