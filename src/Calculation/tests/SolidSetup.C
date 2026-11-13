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
#include <vector>
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
import qchem.ChargeDensity.Seed;   // SeedStrategy (the Auto seed rule)
import qchem.CMath;                // Pi (the k-point spacing)

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

// THE SHELL WHAT-IF measures as the vet stage does (one metric): the trim autoTrim froze, asked of the UNTRIMMED basis's shell table,
// gives the vet stage's own verdict -- its lambda_min AFTER, at the same k, and no AO drop at its orthoTol; with every shell kept, its
// lambda_min BEFORE and drops.  The table's autoTrim reference is the same trim, and its record is filed one level in.
TEST(SolidSetup, TheShellWhatIfMeasuresTheFrozenTrimAsTheVetStageDid)
{
    UnitCell cell=TwoNa(5.0);
    Lattice_3D lat(cell, ivec3_t(2,2,2));
    const report::Run run("shell what-if");
    const auto setup=SolidSetup::Resolve(lat, GaussianBasisRecipe{.data=BasisSetData::VALENCE_LOWQ_VA, .autoTrim=true}, kNa2);
    const auto vet=report::CurrentRunReport().at("basisTrim");   // a COPY: the scopes below grow the run document
    ASSERT_FALSE(setup->recipe->trim.empty()) << "the compressed cell must be near-dependent, or this test tests nothing";
    const auto table=setup->MakeShellConditioning();
    const auto kept=table->KeptBy(setup->recipe->trim);
    BasisSet::Lattice::ShellConditioning::Answer trimmed, untrimmed;
    { const report::Scope s("trimmed");   trimmed  =table->Analyze(kept, kNa2.orthoTol); }
    { const report::Scope s("untrimmed"); untrimmed=table->Analyze(std::vector<bool>(table->Shells().size(), true), kNa2.orthoTol); }
    const double after=vet.at("lambdaMinAfter"), before=vet.at("lambdaMinBefore");
    // eigenvalue roundoff is eps*||S|| (||S|| ~ 1 here), not relative to lambda_min: an ABSOLUTE tolerance
    EXPECT_NEAR(trimmed.lambdaMin, after, 1e-12) << "deleting rows/columns of S IS the trimmed basis's overlap";
    EXPECT_EQ(trimmed.worstK, vet.at("lambdaMinAfterAt").get<std::string>());
    EXPECT_TRUE(trimmed.Healthy()) << "the vet stage stopped where no k-block drops an AO at its orthoTol";
    EXPECT_NEAR(untrimmed.lambdaMin, before, 1e-12);
    EXPECT_FALSE(untrimmed.Healthy());
    EXPECT_EQ(table->KeptBy(setup->MakeVetTrim(kNa2.orthoTol)), kept) << "the reference trim is the loop Resolve ran";
    double share=0.0; for (double p : trimmed.participation) share+=p;
    EXPECT_NEAR(share, 1.0, 1e-10) << "the near-null space is shared out over the kept shells";
    const auto rec=report::CurrentRunReport().at("trimmed");
    EXPECT_EQ(rec.at("shells").size(), table->Shells().size());
    EXPECT_EQ(rec.at("blocks").size(), trimmed.blocks.size());
    EXPECT_TRUE(rec.at("healthy").get<bool>());
}

TEST(SolidSetup, TheShellWhatIfNeedsARecipe)
{
    UnitCell cell=TwoNa(5.0);
    Lattice_3D lat(cell, ivec3_t(1,1,1));
    std::shared_ptr<const BasisSet::Real_BS> mol(BasisSet::Gaussian::Factory(BasisSetData::VALENCE_LOWQ_VA, &cell,
        BasisSet::Gaussian::Engine::MnD, BasisSet::Gaussian::Angular::Cartesian, ShellTrim{}));
    const auto setup=SolidSetup::Resolve(lat, mol, kNa2);
    EXPECT_THROW(setup->MakeShellConditioning(), std::invalid_argument) << "a pre-built basis has no recipe to clear the trim of";
    EXPECT_THROW(setup->MakeVetTrim(1e-4), std::invalid_argument);
}

// The vet stage and the what-if build their Bloch blocks over a basis the recipe made at the RUN's pair-loop tolerances, so they
// must ask for those same tolerances: a stated screenEps once made every autoTrim deck throw GPW_Evaluator's mismatch (2026-11-11).
TEST(SolidSetup, AnAutoTrimAtStatedPairLoopTolerancesResolves)
{
    UnitCell cell=TwoNa(5.0);
    Lattice_3D lat(cell, ivec3_t(2,2,2));
    SolidCalcOptions o=kNa2;
    o.tolerances.screenEps=1e-13;
    const auto setup=SolidSetup::Resolve(lat, GaussianBasisRecipe{.data=BasisSetData::VALENCE_LOWQ_VA, .autoTrim=true}, o);
    EXPECT_FALSE(setup->recipe->trim.empty());
    EXPECT_FALSE(setup->MakeShellConditioning()->Shells().empty());
    EXPECT_EQ(setup->MakeVetTrim(o.orthoTol).shells.size(), setup->recipe->trim.shells.size());
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
enum class XCRhoCell { Becke, GPWRaster, PlaneWaveRaster };   // the solid cells the rho axis splits on
std::vector<std::string> XCRhoRows(XCRhoCell cell, Hamiltonian::XCρ rho, const std::string& tag)
{
    const auto mat=Materials::Get("Si_diamond");
    SolidCalcOptions o; o.Nelec=8; o.densityEcut=20.0; o.accuracy.eps=0.0; o.xcRho=rho;
    if (cell==XCRhoCell::Becke) o.xcMesh=qcMesh::BeckeXCParams(20, -1.0, 11);   // a small recipe: only the setup runs
    else                        o.xcMesh.cellKind=qcMesh::UnitCellKind::Uniform;
    if (cell==XCRhoCell::PlaneWaveRaster) { o.orbitalBasis=SolidCalcOptions::OrbitalBasis::PlaneWave; o.planeWaveEcut=4.0; }
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
    const auto becke=XCRhoRows(XCRhoCell::Becke, Hamiltonian::XCρ::Fit, "Si rho fit Becke");
    ASSERT_EQ(becke.size(), 1u);
    EXPECT_NE(becke[0].find("fit is exact here"), std::string::npos) << becke[0];
    // Plane-wave orbitals: the Vxc ball is >= 4 Ecut, so it holds every orbital pair and P rho = rho[D] -- the row says so.
    const auto pw=XCRhoRows(XCRhoCell::PlaneWaveRaster, Hamiltonian::XCρ::Fit, "Si rho fit PW raster");
    ASSERT_EQ(pw.size(), 1u);
    EXPECT_NE(pw[0].find("fit is exact here"), std::string::npos) << pw[0];
    EXPECT_TRUE(XCRhoRows(XCRhoCell::GPWRaster, Hamiltonian::XCρ::Exact, "Si rho exact raster").empty()) << "the default files no rho row";
}

// GPW orbitals on the plane-wave Vxc fit: no ball holds a Gaussian pair product, and the realized route collocates rho[D] -- so a stated
// xc.rho fit has no realization there.  It must be REFUSED by name, never run as rho[D] (the silent substitution; controller 2026-11-10).
TEST(SolidSetup, AStatedXCRhoFitOnTheGPWBallIsRefusedByName)
{
    try { (void)XCRhoRows(XCRhoCell::GPWRaster, Hamiltonian::XCρ::Fit, "Si rho fit GPW raster"); FAIL() << "xc.rho fit on the GPW ball ran"; }
    catch (const std::invalid_argument& e)
    {
        const std::string what=e.what();
        EXPECT_NE(what.find("xc.rho fit on the plane-wave Vxc fit over Gaussian (GPW) orbitals"), std::string::npos) << what;
        EXPECT_NE(what.find("state xc.rho exact"), std::string::npos) << what;
    }
}

// THE AUTO SEED RULE ON A LATTICE (user rulings 2026-11-12, qchem.SeedPolicy): the default seed is resolved into the setup -- a lattice whose
// species the valence library holds takes IonicSAD (its Madelung field holds the ionic charge states); a species the library lacks takes the
// core guess with a WARN row, where until 2026-11-12 it threw at seed time (Al has no library entry).  A stated seed is honoured as stated.
namespace
{
//! The \a name policy rows a run \a tag filed (anywhere under it).
std::vector<report::json> PolicyRows(const std::string& tag, const std::string& name)
{
    const auto& G=report::GlobalReport();
    const report::json* filed=nullptr;
    for (auto it=G.begin(); it!=G.end(); ++it) if (it.key().rfind(tag+"@",0)==0) filed=&it.value();
    if (!filed) throw std::logic_error(tag+" was not filed");
    std::vector<report::json> out;
    std::function<void(const report::json&)> walk=[&](const report::json& x)
    {
        if (x.is_object())
            for (auto it=x.begin(); it!=x.end(); ++it)
            {
                if (it.key()=="policy" && it.value().is_array())
                    for (const auto& r : it.value()) if (r.value("name","")==name) out.push_back(r);
                walk(it.value());
            }
        else if (x.is_array()) for (const auto& y : x) walk(y);
    };
    walk(*filed);
    return out;
}

std::vector<report::json> AutoSeedRows(const std::string& material, BasisSetData data, ChargeDensity::SeedStrategy seed,
                                       const std::string& tag, ChargeDensity::SeedStrategy& resolved)
{
    const auto mat=Materials::Get(material);
    SolidCalcOptions o; o.densityEcut=20.0; o.accuracy.eps=0.0; o.xcMesh.cellKind=qcMesh::UnitCellKind::Uniform; o.seed=seed;
    { report::Run run(tag);
      resolved=SolidSetup::Resolve(Lattice_3D(*mat.cell, ivec3_t(1,1,1)), GaussianBasisRecipe{.data=data}, o)->opts.seed; }
    return PolicyRows(tag, "AutoSeed");
}
} // namespace

// ONE SOURCE FOR THE MAGNETIC ORDERING (OpenWork §3, 2026-11-12): a stated siteSpins is what the seed reads (the seed key carries
// it), and where its flips differ from the structure's -- an inline cell sent without atom spins -- the setup reports the override.
TEST(SolidSetup, AStatedMagneticOrderingIsTheSeedsAndOverridesTheStructuresSpins)
{
    const auto mno=Materials::Get("MnO_AFM2");
    UnitCell spinless(mno.cell->GetCellMatrix());
    for (Site s : mno.cell->Sites()) { s.spinFlip=false; spinless.Insert(s); }
    auto resolve=[](const UnitCell& cell, std::vector<int> spins, const std::string& tag)
    {
        SolidCalcOptions o; o.densityEcut=20.0; o.accuracy.eps=0.0; o.xcMesh.cellKind=qcMesh::UnitCellKind::Uniform;
        o.multiplicity=1; o.imposeSymmetry=true; o.siteSpins=std::move(spins);
        report::Run run(tag);
        return SolidSetup::Resolve(Lattice_3D(cell, ivec3_t(1,1,1)), GaussianBasisRecipe{.data=BasisSetData::VALENCE_LOWQ_SR}, o)->seedKey.siteSpins;
    };
    const std::vector<int> afm{+1,-1,0,0};
    EXPECT_EQ(resolve(spinless, afm, "MnO spinless stated"), afm) << "the seed key carries the stated ordering";
    const auto rows=PolicyRows("MnO spinless stated", "SeedSpins");
    ASSERT_EQ(rows.size(), 1u) << "a stated ordering whose flips differ from the structure's is reported";
    EXPECT_NE(rows[0].value("decided","").find("override"), std::string::npos) << rows[0].dump();
    EXPECT_NE(rows[0].value("decided","").find("Mn2(-)"), std::string::npos) << "it names the site that differs: " << rows[0].dump();

    EXPECT_EQ(resolve(*mno.cell, afm, "MnO named stated"), afm);
    EXPECT_TRUE(PolicyRows("MnO named stated", "SeedSpins").empty()) << "the material's own spins agree: nothing to report";
}

TEST(SolidSetup, TheAutoSeedIsIonicOnACoveredLatticeAndTheCoreGuessWithAWarnOnAnUncoveredOne)
{
    using ChargeDensity::SeedStrategy;
    SeedStrategy s=SeedStrategy::Default;
    const auto si=AutoSeedRows("Si_diamond", BasisSetData::DZVP, SeedStrategy::Default, "Si auto seed", s);
    EXPECT_EQ(s, SeedStrategy::IonicSAD) << "a covered lattice";
    ASSERT_EQ(si.size(), 1u);
    EXPECT_EQ(si[0].value("verdict",""), "choose");
    EXPECT_NE(si[0].value("reason","").find("Madelung"), std::string::npos) << si[0].dump();

    const auto al=AutoSeedRows("Al_fcc", BasisSetData::VALENCE_LOWQ_SR, SeedStrategy::Default, "Al auto seed", s);
    EXPECT_EQ(s, SeedStrategy::CoreGuess) << "Al has no valence-library entry: the core guess, not a throw";
    ASSERT_EQ(al.size(), 1u);
    EXPECT_EQ(al[0].value("verdict",""), "warn");
    EXPECT_NE(al[0].value("decided","").find("Al"), std::string::npos) << "the warning names the missing species: " << al[0].dump();

    EXPECT_TRUE(AutoSeedRows("Al_fcc", BasisSetData::VALENCE_LOWQ_SR, SeedStrategy::Uniform, "Al stated seed", s).empty()) << "a stated seed files no Auto row";
    EXPECT_EQ(s, SeedStrategy::Uniform);
}

// THE k-SAMPLING POLICY (controller ruling 2026-11-13, from the NiO GUI report): every solid run files ONE `KSampling` row stating the mesh's
// spacing (rigorous: max_i |b_i|/N_i) and the calibrated threshold; a coarser mesh is a WARN naming the mesh that meets the threshold.
// Si's FCC primitive cell: |b_i| = 2π√3/a = 1.0607 bohr^-1 at a = 10.26, so Γ (1.061) and 2x2x2 (0.530) are coarse, 3x3x3 (0.354) is not.
TEST(SolidSetup, AKMeshCoarserThanTheCalibratedSpacingIsWarnedNamingTheMeshThatMeetsIt)
{
    const auto si=Materials::Get("Si_diamond");
    auto rows=[&](const ivec3_t& k, const std::string& tag)
    {
        SolidCalcOptions o; o.densityEcut=20.0; o.accuracy.eps=0.0; o.xcMesh.cellKind=qcMesh::UnitCellKind::Uniform;
        { report::Run run(tag); (void)SolidSetup::Resolve(Lattice_3D(*si.cell, k), GaussianBasisRecipe{.data=BasisSetData::SIPP}, o); }
        return PolicyRows(tag, "KSampling");
    };
    const auto g=rows(ivec3_t(1,1,1), "Si k111");
    ASSERT_EQ(g.size(), 1u) << "one KSampling row per run";
    EXPECT_EQ(g[0].value("verdict",""), "warn");
    const double b=2*Pi*si.cell->InverseInterplanarSpacings().x;
    EXPECT_NEAR(g[0]["detail"]["spacing"].get<double>(), b, 1e-12) << "Γ only: the spacing is |b| itself";
    EXPECT_GT(g[0]["detail"]["threshold"].get<double>(), 0.0);
    EXPECT_NE(g[0].value("decided","").find("a 3x3x3 mesh meets it"), std::string::npos) << g[0].dump();

    const auto k2=rows(ivec3_t(2,2,2), "Si k222");
    ASSERT_EQ(k2.size(), 1u);
    EXPECT_EQ(k2[0].value("verdict",""), "warn") << k2[0].dump();
    EXPECT_NEAR(k2[0]["detail"]["spacing"].get<double>(), b/2, 1e-12);

    const auto k3=rows(ivec3_t(3,3,3), "Si k333");
    ASSERT_EQ(k3.size(), 1u);
    EXPECT_EQ(k3[0].value("verdict",""), "choose") << "a mesh within the threshold is reported, not warned: " << k3[0].dump();
}
