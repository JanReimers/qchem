// File: Calculation/tests/SolidSetup.C  The periodic RUN SETUP (doc/CalcNetworkReview.md §7, step 6): its basis recipe
// freezes the vet-stage trim at the reference geometry, and checks the basis's pseudopotential provenance with no deck.
#include <gtest/gtest.h>
#include <unistd.h>    // getpid: per-PROCESS temp names (ctest runs each case as its own process, in parallel)
#include <filesystem>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
import qchem.SolidCalculation;
import qchem.Lattice_3D;
import qchem.Structure;
import qchem.UnitCell;
import qchem.BasisSet;
import qchem.BasisSet.Gaussian.Point.Factory;
import qchem.BasisSet.Lattice.BasisSet;   // VetStageTrim, ShellTrim
import qchem.Types;

using namespace qchem;
using BasisSet::Gaussian::BasisSetData;
using BasisSet::Gaussian::ShellTrim;

namespace
{
//! Two Na in a cubic cell of edge \a a: compressed (a = 5) their diffuse s/p are near-dependent, so the vet loop trims.
UnitCell TwoNa(double a)
{
    UnitCell cell(a);
    cell.AddAtom(11,{0.0,0.0,0.0});
    cell.AddAtom(11,{0.5,0.5,0.5});
    return cell;
}
const SolidCalcOptions kNa2{.Nelec=2, .pp="gth-lda-valence", .orthoTol=1e-4};
} // namespace

TEST(SolidSetup, ARecipeFreezesTheTrimTheVetLoopDecides)
{
    UnitCell cell=TwoNa(5.0);
    Lattice_3D lat(cell, ivec3_t(2,2,2));
    const auto setup=SolidSetup::Resolve(lat, GaussianBasisRecipe{.data=BasisSetData::VALENCE_LOWQ_VA, .vet=true}, kNa2);
    ASSERT_TRUE(setup->recipe);
    EXPECT_FALSE(setup->recipe->vet) << "the vet loop runs ONCE, at the reference geometry; the recipe keeps its decision";
    // The deck's twin: the same loop, run by hand on the same cell.
    auto make=[&](const ShellTrim& t){ return std::shared_ptr<const BasisSet::Real_BS>(BasisSet::Gaussian::Factory(
        BasisSetData::VALENCE_LOWQ_VA, &cell, BasisSet::Gaussian::Engine::MnD, BasisSet::Gaussian::Angular::Cartesian, t)); };
    const auto twin=BasisSet::Lattice::VetStageTrim(lat, make, {}, kNa2.orthoTol);
    ASSERT_FALSE(twin.trim.empty()) << "the compressed cell must be near-dependent, or this test tests nothing";
    const ShellTrim& frozen=setup->recipe->trim;
    ASSERT_EQ(frozen.shells.size(), twin.trim.shells.size());
    for (const auto& s : twin.trim.shells)
        EXPECT_TRUE(frozen.Removes(s.Z, s.l, s.exponents)) << "Z=" << s.Z << " l=" << s.l;
    EXPECT_EQ(setup->mol->GetNumFunctions(), twin.mol->GetNumFunctions());
    // Rebuilt from the frozen recipe -- what a geometry step does -- the basis is the trimmed one, with no second vet.
    EXPECT_EQ(setup->recipe->Make(*setup->cell, setup->pp)->GetNumFunctions(), twin.mol->GetNumFunctions());
}

TEST(SolidSetup, APreBuiltBasisHasNoRecipe)
{
    UnitCell cell=TwoNa(20.0);
    Lattice_3D lat(cell, ivec3_t(1,1,1));
    std::shared_ptr<const BasisSet::Real_BS> mol(BasisSet::Gaussian::Factory(BasisSetData::VALENCE_LOWQ_VA, &cell));
    const auto setup=SolidSetup::Resolve(lat, mol, kNa2);
    EXPECT_FALSE(setup->recipe) << "a caller-built basis cannot be rebuilt at moved atoms";
    EXPECT_EQ(setup->mol, mol);
}

TEST(SolidSetup, ARecipeChecksTheBasisProvenanceWithNoDeck)
{
    const std::string path=(std::filesystem::path(testing::TempDir())/(std::to_string(::getpid())+"_na9.bsd")).string();
    { std::ofstream f(path); f<<"title\nPP-PROVENANCE NA=9\n\n!\nBASIS=\"t\"\n"; }
    UnitCell cell=TwoNa(20.0);
    Lattice_3D lat(cell, ivec3_t(1,1,1));
    try { (void)SolidSetup::Resolve(lat, GaussianBasisRecipe{.file=path}, kNa2); FAIL() << "a q1 run on a q9 file must throw"; }
    catch (const std::runtime_error& e)
    {
        const std::string m=e.what();
        EXPECT_NE(m.find("Na"),std::string::npos) << m;
        EXPECT_NE(m.find("q9"),std::string::npos) << m;
        EXPECT_NE(m.find("q1"),std::string::npos) << m;
    }
}
