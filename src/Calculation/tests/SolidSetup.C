// File: Calculation/tests/SolidSetup.C  The periodic RUN SETUP (doc/CalcNetworkReview.md §7, step 6): its basis recipe
// freezes the vet-stage trim at the reference geometry, and checks the basis's pseudopotential provenance with no deck.
#include <gtest/gtest.h>
#include <unistd.h>    // getpid: per-PROCESS temp names (ctest runs each case as its own process, in parallel)
#include <filesystem>
#include <fstream>
#include <functional>
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
import qchem.BasisSet.Gaussian.Lattice.LatticeSum1E;   // GaussianSharpness::Tolerances (what the built basis computes at)
import qchem.Materials;   // Si_diamond (the xc.rho policy case)
import qchem.Mesh;        // UnitCellKind
import qchem.Mesh.XCPolicy;   // BeckeXCParams
import qchem.Reporting;   // the policy rows the setup files

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
    const auto setup=SolidSetup::Resolve(lat, GaussianBasisRecipe{.data=BasisSetData::VALENCE_LOWQ_VA, .autoTrim=true}, kNa2);
    ASSERT_TRUE(setup->recipe);
    EXPECT_FALSE(setup->recipe->autoTrim) << "the vet loop runs ONCE, at the reference geometry; the recipe keeps its decision";
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

TEST(SolidSetup, AnAutoTrimWithNoOrthoTolIsRefused)
{
    Lattice_3D lat(TwoNa(5.0), ivec3_t(2,2,2));
    SolidCalcOptions o=kNa2; o.orthoTol=0.0;
    EXPECT_THROW(SolidSetup::Resolve(lat, GaussianBasisRecipe{.data=BasisSetData::VALENCE_LOWQ_VA, .autoTrim=true}, o), std::invalid_argument)
        << "at orthoTol 0 it would drop nothing, silently";
}

TEST(SolidSetup, TheRungGateCalibrationReachesTheLadderTolerances)
{
    Lattice_3D lat(TwoNa(5.0), ivec3_t(1,1,1));
    SolidCalcOptions o=kNa2;
    const GaussianBasisRecipe recipe{.data=BasisSetData::VALENCE_LOWQ_VA};
    EXPECT_EQ(SolidSetup::Resolve(lat, recipe, o)->opts.tolerances.rungGateC, BasisSet::Gaussian::kRungGateC);
    o.accuracy.rungGateC=6.0;
    EXPECT_EQ(SolidSetup::Resolve(lat, recipe, o)->opts.tolerances.rungGateC, 6.0);
}

// The pair rule's stiffness (V-BS-RASTER step 2) is ONE value read twice: the ladder's tolerances and the basis's own pair->level
// assignment.  A stated calibration must reach BOTH -- the basis is built before Decide derives the run's tolerances, and a basis bound
// with the default would make the GPW evaluator refuse the run (it checks the two agree).
TEST(SolidSetup, TheRelCutoffSafetyCalibrationReachesTheLadderAndTheBasis)
{
    Lattice_3D lat(TwoNa(5.0), ivec3_t(1,1,1));
    SolidCalcOptions o=kNa2;
    const GaussianBasisRecipe recipe{.data=BasisSetData::VALENCE_LOWQ_VA};
    auto bound=[](const SolidSetup& s)
    {
        for (auto ibs : s.mol->Iterate<BasisSet::Real_OIBS>())
            if (const auto* g=dynamic_cast<const BasisSet::Gaussian::GaussianSharpness*>(ibs)) return g->Tolerances().relCutoffSafety;
        return -1.0;
    };
    const auto d=SolidSetup::Resolve(lat, recipe, o);
    EXPECT_EQ(d->opts.tolerances.relCutoffSafety, BasisSet::Gaussian::kRelCutoffSafety);
    EXPECT_EQ(bound(*d), BasisSet::Gaussian::kRelCutoffSafety);
    o.accuracy.relCutoffSafety=3.0;
    const auto s=SolidSetup::Resolve(lat, recipe, o);
    EXPECT_EQ(s->opts.tolerances.relCutoffSafety, 3.0) << "the ladder's";
    EXPECT_EQ(bound(*s), 3.0) << "the basis's pair rule";
    o.accuracy.relCutoffSafety=0.0;
    EXPECT_THROW(SolidSetup::Resolve(lat, recipe, o), std::invalid_argument) << "a stiffness is > 0";
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

// D-ENV-TOL option A (2026-10-31): the recipe BINDS the run's pair-loop tolerances on the basis it builds -- at the reference and
// on every rebuild from the frozen recipe -- and the GPW evaluator only checks them, so a pre-built (default) basis under a run
// that asks for others is REFUSED rather than silently re-tuned.
TEST(SolidSetup, ARecipeBindsTheRunsGPWTolerancesAndAPreBuiltBasisCannotTakeOthers)
{
    UnitCell cell=TwoNa(20.0);
    Lattice_3D lat(cell, ivec3_t(1,1,1));
    SolidCalcOptions o=kNa2; o.tolerances.screenEps=1e-8; o.tolerances.densityEps=1e-11;
    auto bound=[](const BasisSet::Real_BS& bs)
    {
        for (auto ibs : bs.Iterate<BasisSet::Real_OIBS>())
            if (auto* g=dynamic_cast<const BasisSet::Gaussian::GaussianSharpness*>(ibs)) return g->Tolerances();
        ADD_FAILURE() << "no Gaussian block"; return BasisSet::Gaussian::GPWTolerances{};
    };
    const auto setup=SolidSetup::Resolve(lat, GaussianBasisRecipe{.data=BasisSetData::VALENCE_LOWQ_VA}, o);
    ASSERT_TRUE(setup->recipe);
    EXPECT_EQ(setup->recipe->tolerances.screenEps, 1e-8) << "frozen with the trim";
    EXPECT_EQ(bound(*setup->mol).screenEps, 1e-8);
    EXPECT_EQ(bound(*setup->mol).densityEps, 1e-11);
    EXPECT_EQ(bound(*setup->recipe->Make(*setup->cell, setup->pp)).screenEps, 1e-8) << "a geometry step's basis computes at the same setting";

    std::shared_ptr<const BasisSet::Real_BS> mol(BasisSet::Gaussian::Factory(BasisSetData::VALENCE_LOWQ_VA, &cell));
    EXPECT_EQ(bound(*mol).screenEps, BasisSet::Gaussian::GPWTolerances{}.screenEps) << "a pre-built basis carries the defaults";
    EXPECT_THROW(SolidCalculation(SolidSetup::Resolve(lat, mol, o)), std::invalid_argument)
        << "the run asks for tolerances its pre-built basis was not built with";
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

// THE ELECTRON COUNT IS THE SPECIES' (2026-10-31): a stated count that disagrees with sum(q) - charge is REFUSED, naming both.  The
// old silent +8: a stated SEMICORE Na q9 with the material's q1 count (8 for NaF instead of 16) ran a charged cell without a word.
TEST(SolidSetup, AnElectronCountThatDisagreesWithTheSpeciesIsRefused)
{
    UnitCell cell(8.0);
    cell.AddAtom(11,{0.0,0.0,0.0});
    cell.AddAtom( 9,{0.5,0.5,0.5});
    Lattice_3D lat(cell, ivec3_t(1,1,1));
    const GaussianBasisRecipe q9{.file="valence_lowq_sr2_naq9.bsd"};
    const SolidCalcOptions semicore{.species={{"Na",9},{"F",7}}, .pp="gth-lda-valence"};
    SolidCalcOptions old=semicore; old.Nelec=8;   // the q1 count, the bug
    try { (void)SolidSetup::Resolve(lat, q9, old); FAIL() << "a stated Nelec 8 against Na q9 + F q7 must be refused"; }
    catch (const std::invalid_argument& e)
    {
        const std::string w=e.what();
        EXPECT_NE(w.find("Nelec 8"), std::string::npos) << w;
        EXPECT_NE(w.find("16"), std::string::npos) << "names the species' count: " << w;
    }
    EXPECT_EQ(SolidSetup::Resolve(lat, q9, semicore)->opts.Nelec, 16) << "unstated: derived from the species";
    SolidCalcOptions stated=semicore; stated.Nelec=16;
    EXPECT_NO_THROW((void)SolidSetup::Resolve(lat, q9, stated)) << "a stated count that agrees";
    SolidCalcOptions cation=semicore; cation.charge=2;
    EXPECT_EQ(SolidSetup::Resolve(lat, q9, cation)->opts.Nelec, 14) << "a STATED charge is a legal charged cell";
    cation.Nelec=14;
    EXPECT_NO_THROW((void)SolidSetup::Resolve(lat, q9, cation));
}

namespace
{
//! Resolve Si (Gamma, stated grids) on a Becke mesh or the raster at \a rho in a run named \a tag; the `decided` lines of its XCRhoProjected rows.
std::vector<std::string> XCRhoRows(bool becke, Hamiltonian::XCρ rho, const std::string& tag)
{
    const auto mat=Materials::Get("Si_diamond");
    SolidCalcOptions o; o.Nelec=8; o.densityEcut=20.0; o.accuracy.eps=0.0; o.xcRho=rho;
    if (becke) o.xcMesh=qcMesh::BeckeXCParams(20, -1.0, 11);   // a small recipe: only the setup runs
    else         o.xcMesh.cellKind=qcMesh::UnitCellKind::Uniform;
    { report::Run run(tag);
      (void)SolidSetup::Resolve(Lattice_3D(*mat.cell, ivec3_t(1,1,1)), GaussianBasisRecipe{.data=BasisSetData::SIPP_SR}, o); }
    const auto& G=report::GlobalReport();
    const report::json* filed=nullptr;
    for (auto it=G.begin(); it!=G.end(); ++it) if (it.key().rfind(tag+"@",0)==0) filed=&it.value();
    std::vector<std::string> out;
    if (!filed) throw std::logic_error(tag+" was not filed");
    std::function<void(const report::json&)> walk=[&](const report::json& x)
    {
        if (x.is_object())
            for (auto it=x.begin(); it!=x.end(); ++it)
            {
                if (it.key()=="policy" && it.value().is_array())
                    for (const auto& r : it.value()) if (r.value("name","")=="XCRhoProjected") out.push_back(r.value("decided",""));
                walk(it.value());
            }
        else if (x.is_array()) for (const auto& y : x) walk(y);
    };
    walk(*filed);
    return out;
}
} // namespace

//! THE XC rho AXIS ON A CELL (user ruling 2026-11-09, V-XC-ONE-TERM): a stated `xc.rho fit` is ACCEPTED, and the setup files the
//! choice on the fit the Hamiltonian RESOLVES -- Auto on a Becke mesh is the delta fit ("fit is exact here"), on the raster the
//! plane-wave fit (its G-ball's projection).  The default (exact) files no rho row.
TEST(SolidSetup, AStatedXCRhoFitIsAcceptedAndReportedOnTheResolvedFit)
{
    const auto becke=XCRhoRows(true, Hamiltonian::XCρ::Fit, "Si rho fit Becke");
    ASSERT_EQ(becke.size(), 1u);
    EXPECT_NE(becke[0].find("fit is exact here"), std::string::npos) << becke[0];
    const auto raster=XCRhoRows(false, Hamiltonian::XCρ::Fit, "Si rho fit raster");
    ASSERT_EQ(raster.size(), 1u);
    EXPECT_NE(raster[0].find("plane-wave"), std::string::npos) << raster[0];
    EXPECT_TRUE(XCRhoRows(false, Hamiltonian::XCρ::Exact, "Si rho exact raster").empty()) << "the default files no rho row";
}
