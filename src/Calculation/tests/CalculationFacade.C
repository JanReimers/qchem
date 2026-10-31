// File: Calculation/tests/CalculationFacade.C  What the molecular facade hands a front end (lane G, 2026-10-30): the
// constructor's own SCF observed from its first iteration, every MO as a level, and the spin density -- a zero field when
// the run is unpolarized, never an absence.
#include <gtest/gtest.h>
#include <cmath>
#include <vector>
import qchem.Calculation;
import qchem.StructureData;
import qchem.SCFParams;
import qchem.SCFIterator;
import qchem.Types;

using namespace qchem;
using StructureData::GetMolecule;

TEST(CalculationFacade, TheConstructorsObserverSeesItsSCFFromTheFirstIteration)
{
    std::vector<SCFIterator::SCFProgress> steps;
    Calculation calc(MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp"}), {},
                     [&steps](const SCFIterator::SCFProgress& p){ steps.push_back(p); });
    ASSERT_TRUE(calc.IsConverged());
    ASSERT_FALSE(steps.empty()) << "the observer must see the constructor's SCF, not only a later Converge";
    EXPECT_EQ(steps.front().iteration, 1u) << "iteration 1: the constructor's SCF from its start";
    EXPECT_EQ(steps.size(), calc.IterationCount()) << "every iteration, none skipped";
    EXPECT_NEAR(steps.back().energy, calc.Energy(), 1e-10);
}

TEST(CalculationFacade, LevelsListEveryMOAscendingAndTheOccupiedOnesAreTheOrbitals)
{
    Calculation calc(MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp"}));
    ASSERT_TRUE(calc.IsConverged());
    const auto& levels=calc.GetLevels();
    ASSERT_GT(levels.size(), calc.NumOccupied()) << "the virtuals are listed too";
    double ne=0; size_t nocc=0;
    for (size_t i=0;i<levels.size();i++)
    {
        ASSERT_NE(levels[i].phi, nullptr);
        if (i>0) EXPECT_LE(levels[i-1].e, levels[i].e) << "ascending";
        ne+=levels[i].occupation;
        if (levels[i].occupation>0.5) nocc++;
        EXPECT_EQ(levels[i].spin, Spin::None) << "an unpolarized run";
    }
    EXPECT_NEAR(ne, 10.0, 1e-12);
    ASSERT_EQ(nocc, calc.NumOccupied());
    EXPECT_EQ(levels[nocc-1].phi, &calc.HOMO()) << "the highest occupied level IS the facade's HOMO";
    EXPECT_GT(levels[nocc].e, levels[nocc-1].e) << "a gap to the LUMO";
}

TEST(CalculationFacade, AnUnpolarizedSpinDensityIsTheZeroFieldAndAPolarizedOneIsNot)
{
    Calculation closed(MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp"}));
    ASSERT_TRUE(closed.IsConverged());
    EXPECT_EQ(closed.GetSpinDensity()(rvec3_t(0.0, 0.0, 0.3)), 0.0);

    Calculation triplet(MoleculeSetup::Resolve(GetMolecule("O2"), {.basis="dzvp", .multiplicity=3}));
    ASSERT_TRUE(triplet.IsConverged());
    // The two unpaired π* electrons put the moment off the bond axis, around the atoms.
    EXPECT_GT(triplet.GetSpinDensity()(rvec3_t(0.8, 0.0, 1.141)), 1e-3);
    size_t up=0, down=0; double nUp=0, nDown=0;
    for (const auto& L : triplet.GetLevels())
    {
        if (L.spin==Spin::Up)   { up++;   nUp+=L.occupation; }
        if (L.spin==Spin::Down) { down++; nDown+=L.occupation; }
    }
    EXPECT_EQ(up, down) << "both channels listed";
    EXPECT_NEAR(nUp-nDown, 2.0, 1e-12) << "2S = multiplicity - 1";
}
