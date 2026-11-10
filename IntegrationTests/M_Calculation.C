// File: UnitTests/M_Calculation.C  Tests for the qchem::MoleculeCalculation front-door facade.
//
// The facade is meant to BE the production recipe -- so its job here is to reproduce, through the
// public front door alone, the converged HF/dzvp water energy that M_HF_U checks via the QchemTester
// scaffold.  Same geometry, same basis, same anchor.  Also a quick smoke of the ScalarFunction
// sampling surface (Density/HOMO/Orbital), the win the API-ergonomics review highlighted.
#include "gtest/gtest.h"
#include <unistd.h>    // getpid: per-process temp names (ctest runs cases in parallel)
#include <cmath>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

import qchem.MoleculeCalculation;
import qchem.AtomCalculation;     // a molecule state offered to an atom run (refused by kind)
import qchem.SCFIterator;         // IterationRow (a restart's iterates)
import qchem.SCFParams;
import qchem.ChargeDensity.Seed;   // SeedStrategy (HF from a SAD seed)
import qchem.StructureData;        // GetMolecule("H2O") -- the shared water geometry (D-MAKEWATER)
import qchem.Structure;
import qchem.Types;        // Vector3D
import qchem.Reporting;    // report:: -- the run report the facade emits (scf section)
using namespace qchem;

using qchem::MoleculeCalculation;

// Regression anchor: identical converged HF/dzvp total energy to M_HF_U_Water (-76.022903).
TEST(M_Calculation, WaterEnergy)
{
    Molecule water = qchem::StructureData::GetMolecule("H2O");
    MoleculeCalculation calc(MoleculeSetup::Resolve(water, {.basis = "dzvp"}));            // build + converge in one line

    const double E_ref = -76.022903;
    EXPECT_LT(std::fabs((E_ref - calc.Energy()) / E_ref), 1e-5);
    EXPECT_GT(calc.IterationCount(), 0u);
}

// The scf run-report SCHEMA CHECK (RunReportPlan step 1).  The provider (MoleculeCalculation::Converge)
// builds the `scf` section json directly; here we drive one Converge and assert the emitted shape --
// this is the typo-catcher that replaces the compile-time key-safety we traded for "the report IS json".
TEST(M_Calculation, ScfReportSchema)
{
    Molecule water = qchem::StructureData::GetMolecule("H2O");
    { report::Run run(water.ID()); MoleculeCalculation calc(MoleculeSetup::Resolve(water, {.basis = "dzvp"})); }   // the test opens the run (user ruling D, 2026-11-02); the ctor's Converge files scf into it

    // Find this run's scf section (the run key carries a timestamp, so scan for it).
    const report::json* scf = nullptr;
    const report::json& all = report::GlobalReport();
    for (auto it = all.begin(); it != all.end(); ++it)
        if (it.value().contains("scf")) scf = &it.value()["scf"];
    ASSERT_NE(scf, nullptr) << "no run emitted an scf section";

    ASSERT_TRUE(scf->contains("standard"));
    ASSERT_TRUE(scf->contains("advanced"));
    const report::json& std_ = (*scf)["standard"];
    for (const char* k : { "nMaxIter", "minDrho", "minDFD", "minDE", "minFD",
                           "mixer", "kerkerG0", "accelerator", "smearingkT" })
        EXPECT_TRUE(std_.contains(k)) << "scf.standard missing key: " << k;
    EXPECT_EQ(std_["accelerator"].get<std::string>(), "DIIS");   // facade default
    EXPECT_EQ(std_["mixer"].get<std::string>(), "Lin");          // no Kerker/Pulay by default

    const report::json& adv = (*scf)["advanced"];
    for (const char* k : { "pulayDepth", "pulayStart", "useMOM", "momStartIter", "momSmearPenalty" })
        EXPECT_TRUE(adv.contains(k)) << "scf.advanced missing key: " << k;
    ASSERT_TRUE(adv.contains("guard"));
    EXPECT_TRUE(adv["guard"].contains("holePersistence"));
    EXPECT_TRUE(adv["guard"].contains("maxReleases"));
}

// The basis run-report SCHEMA CHECK (RunReportPlan step 2).  This is the cross-layer / cursor path:
// the orchestrator stamps the scalars, while basis.perIrrep is filled by the LASolver's conditioning
// write five layers down (MakeIrrepWFs owns the row).  Symmetry ON so C2v water yields several irreps
// -> perIrrep is a genuine multi-row table, and each row must carry the LASolver-supplied conditioning.
TEST(M_Calculation, BasisReportSchema)
{
    report::ClearGlobal();                                // isolate: only this test's run(s) remain
    { report::Run run("H2O"); MoleculeCalculation calc(MoleculeSetup::Resolve(qchem::StructureData::GetMolecule("H2O"), { .basis = "dzvp", .symmetry = true })); }

    const report::json& all = report::GlobalReport();
    const report::json* basis = nullptr;
    for (auto it = all.begin(); it != all.end(); ++it)
        if (it.value().contains("basis")) basis = &it.value()["basis"];
    ASSERT_NE(basis, nullptr) << "no run emitted a basis section";

    EXPECT_EQ((*basis)["name"].get<std::string>(), "dzvp");
    EXPECT_EQ((*basis)["engine"].get<std::string>(), "mnd");
    EXPECT_EQ((*basis)["angular"].get<std::string>(), "cartesian");
    EXPECT_GT((*basis)["nFunctions"].get<int>(), 0);

    ASSERT_TRUE(basis->contains("perIrrep"));
    ASSERT_TRUE((*basis)["perIrrep"].is_array());
    EXPECT_GE((*basis)["perIrrep"].size(), 2u) << "C2v water should split into several irreps";
    for (const report::json& row : (*basis)["perIrrep"])
        for (const char* k : { "irrep", "nFunctions", "lambdaMin", "lambdaMax", "cond" })
            EXPECT_TRUE(row.contains(k)) << "perIrrep row missing key: " << k;   // cond came from the LASolver
}

// The cache run-report SCHEMA CHECK (RunReportPlan step 4).  Converge snapshots the integrals cache into the
// run's `cache` section AFTER Iterate (inside the run bracket).  The cache is a process-wide singleton never
// cleared between runs, so the numbers are cumulative -- here we only assert the shape.
TEST(M_Calculation, CacheReportSchema)
{
    report::ClearGlobal();                               // isolate the REPORT (the cache itself is never cleared)
    { report::Run run("H2O"); MoleculeCalculation calc(MoleculeSetup::Resolve(qchem::StructureData::GetMolecule("H2O"), { .basis = "dzvp" })); }  // HF over dzvp populates Jac/Kab + the Omega/Hermite caches

    const report::json& all = report::GlobalReport();
    const report::json* cache = nullptr;
    for (auto it = all.begin(); it != all.end(); ++it)
        if (it.value().contains("cache")) cache = &it.value()["cache"];
    ASSERT_NE(cache, nullptr) << "no run emitted a cache section";

    ASSERT_TRUE(cache->contains("tiers"));
    ASSERT_TRUE((*cache)["tiers"].is_array());
    EXPECT_GE((*cache)["tiers"].size(), 2u);            // Jac, Kab (Cach4 unused by HF -> dropped, not shown)
    for (const report::json& t : (*cache)["tiers"])
        for (const char* k : { "name", "ramMB", "pct" })
            EXPECT_TRUE(t.contains(k)) << "tier row missing key: " << k;

    ASSERT_TRUE(cache->contains("reuse"));              // HF/dzvp exercises the Omega + Hermite caches
    ASSERT_TRUE((*cache)["reuse"].is_array());
    ASSERT_GE((*cache)["reuse"].size(), 1u);
    for (const char* k : { "cache", "entries", "lookups", "reusePct", "ramMB" })
        EXPECT_TRUE((*cache)["reuse"][0].contains(k)) << "reuse row missing key: " << k;
}

// Parameter-free LDA molecular DFT through the facade (Dirac exchange + VWN5 correlation).  The facade
// auto-selects the SAD seed and DIIS-from-start for DFT; CalcOptions.mesh defaults to molecular values.
// "Did E move" regression sentinel, like the HF anchor (NOT a physical-accuracy claim).
TEST(M_Calculation, WaterLDA)
{
    MoleculeCalculation calc(MoleculeSetup::Resolve(qchem::StructureData::GetMolecule("H2O"), {.basis = "dzvp", .model = qchem::Model::LDA}));
    // Converged parameter-free LDA total energy (Dirac exchange + VWN5 correlation), with the facade's
    // auto DIIS-from-start + SAD seed.  Bit-stable at -75.9324615507 and now ORDER-INDEPENDENT: the
    // ~585 ppm HF-before-DFT drift was a fit-basis Normalization cached without a mesh key (the HF SAD
    // bootstrap's coarse seed mesh poisoned the production run); fixed by mesh-keying Fit_IBS::Norm().
    // So the anchor can stay tight (1e-5) -- if it moves, something regressed.
    // THE LDA MOLECULE ANCHOR AT A CONVERGED MESH (2026-10-31): the default mesh is MHL 80 / Lebedev 29, the default Vxc fit the
    // δ basis (user ruling "δ default, A1 kept"), and the fit bases keep their D shells.  Was -75.9324615507 on the 30/5 mesh
    // (60 mHa LOW) and -75.8726517 on the A1 Vxc fit with a d-less Coulomb fit.  Independent route: PySCF with the SAME A1_coul
    // Coulomb fit and exact XC, -75.879447 (8 μHa); PySCF exact J -75.877299 (the Coulomb fit itself lowers 2.1 mHa).
    const double E_ref = -75.8794551388;
    EXPECT_LT(std::fabs((E_ref - calc.Energy()) / E_ref), 1e-5);
    EXPECT_NEAR(calc.Energy(), -75.87944736, 1e-4) << "ORACLE: PySCF lda,vwn/dzvp, the same A1_coul Coulomb fit, exact XC";
    EXPECT_GT(calc.IterationCount(), 0u);
}

// {.symmetry=true} SALC-blocks the (Cartesian PG) basis + does global aufbau across irreps.  The
// converged total energy must match the un-blocked run -- the M_Sym invariant, now through the facade.
TEST(M_Calculation, WaterSymmetry)
{
    MoleculeCalculation calc(MoleculeSetup::Resolve(qchem::StructureData::GetMolecule("H2O"), {.basis = "dzvp", .symmetry = true}));

    const double E_ref = -76.022903;                       // identical to the un-blocked WaterEnergy anchor
    EXPECT_LT(std::fabs((E_ref - calc.Energy()) / E_ref), 1e-5);
    EXPECT_GT(calc.IterationCount(), 0u);
}

// engine=LibCint + {.symmetry=true}: the libcint-Cartesian orbital IBS is ALSO symmetry-adaptable (it
// implements GetAoShells via Gaussian::Orbital_1E_IBS, like the MnD bases), so this must converge to the
// same energy as the MnD/un-blocked run.  Guards the {engine=libcint}x{symmetry} combination, previously
// untested -- a symmetry-adaptation refactor once silently regressed it (libcint adapted via the old PGData
// cast); this test locks it in.
TEST(M_Calculation, WaterSymmetryLibCint)
{
    MoleculeCalculation calc(MoleculeSetup::Resolve(qchem::StructureData::GetMolecule("H2O"), {.basis = "dzvp", .engine = Engine::LibCint, .symmetry = true}));

    const double E_ref = -76.022903;                       // same anchor as WaterSymmetry (MnD)
    EXPECT_LT(std::fabs((E_ref - calc.Energy()) / E_ref), 1e-5);
    EXPECT_GT(calc.IterationCount(), 0u);
}

// The caller's Molecule is deep-copied: it survives being passed in and may be used afterwards.
TEST(M_Calculation, OwnsItsOwnStructure)
{
    MoleculeCalculation calc(MoleculeSetup::Resolve(qchem::StructureData::GetMolecule("H2O"), {.basis = "dzvp"}));      // temporary destroyed at the semicolon
    EXPECT_EQ(calc.GetStructure().GetNumAtoms(), 3u);
    EXPECT_NEAR(calc.GetStructure().GetNumElectrons(), 10.0, 1e-12);
}

// Density and every occupied MO are sampleable ScalarFunction<double>s through one interface.
TEST(M_Calculation, SamplingSurface)
{
    MoleculeCalculation calc(MoleculeSetup::Resolve(qchem::StructureData::GetMolecule("H2O"), {.basis = "dzvp"}));

    // Water has 10 electrons -> 5 doubly-occupied MOs in the unpolarized wave function.
    EXPECT_EQ(calc.NumOccupied(), 5u);

    const Vector3D<double> rO(0, 0, 0);                    // at the oxygen nucleus
    EXPECT_GT(calc.Density()(rO), 0.0);                    // rho > 0 where the atoms are
    EXPECT_GT(std::fabs(calc.HOMO()(rO)) + std::fabs(calc.Orbital(0)(rO)), 0.0);
}

// ===== The GDM/Ladder molecular path, stressed with BORON (feedback pin: B, not Be) =================
// Boron's open 2p^1 doublet is the light-atom SCF stress case: plain aufbau+DIIS limit-cycles on the
// degenerate p shell, which is exactly the regime the DIIS->GDM ladder exists for -- and until this
// test the molecular GDM/Ladder path had NO gtest coverage at all (it was reachable only through the
// facade's accelerator option).  Runs UHF (multiplicity=2 promotes to Polarized) on dzvp through the
// public front door with {.type="Ladder"}, and pins the converged energy as a did-E-move anchor.
// This is ALSO the behavioural anchor for the scalar-agnostic accelerator-manager face
// (doc/RealComplexPlan.md §6): one non-template manager driving the real molecular path.
static Molecule MakeBoron()
{
    Molecule b;
    b.Insert(Site{5, Vector3D<double>(0, 0, 0)});
    return b;
}

TEST(M_Calculation, BoronUHFLadderGDM)
{
    MoleculeCalculation calc(MoleculeSetup::Resolve(MakeBoron(), {.basis = "dzvp", .multiplicity = 2}, {.type = "Ladder"}));

    const double E_ref = -24.525932;   // pinned converged UHF/dzvp did-E-move anchor (measured 2026-08-17;
                                       // ~3 mHa above the HF limit -24.529061 = dzvp basis incompleteness)
    EXPECT_LT(std::fabs((E_ref - calc.Energy()) / std::fabs(E_ref)), 1e-5);
    EXPECT_GT(calc.IterationCount(), 0u);
}

// The same anchor under standalone GDM (no DIIS rung in front): the direct minimizer alone must land
// on the same UHF energy.  Locks the second previously-uncovered accelerator type.
TEST(M_Calculation, BoronUHFPureGDM)
{
    MoleculeCalculation calc(MoleculeSetup::Resolve(MakeBoron(), {.basis = "dzvp", .multiplicity = 2}, {.type = "GDM"}));

    const double E_ref = -24.525932;   // same anchor as the Ladder run (the two agree to ~1e-12 Ha)
    EXPECT_LT(std::fabs((E_ref - calc.Energy()) / std::fabs(E_ref)), 1e-5);
    EXPECT_GT(calc.IterationCount(), 0u);
}

// HF with degenerateLevelFill: Shared STATED, from a SAD seed: the symmetric (shared) UHF state, -24.40769 on dzvp, 118 mHa above the
// integer state (the same-spin self-interaction a fractional HF occupation pays).  A matrix-free (SAD) seed reaches HF through an LDA
// sibling whose own fill is integer; the stated rule must refill that seed, or the run starts symmetry-broken and STAYS there.
// DISCRIMINATES (lane A 2026-11-09): without the sibling's refill this run lands on the integer state -24.5259316 (8 iterations).
// The seed must be stated: HF's default seed is the core guess, which never takes the sibling path.
TEST(M_Calculation, BoronUHFStatedSharedTakesTheSharedState)
{
    SCFParams shared; shared.DegenerateLevels=SCFParams::DegenerateLevelFill::Shared;
    MoleculeCalculation calc(MoleculeSetup::Resolve(MakeBoron(), {.basis = "dzvp", .multiplicity = 2,
                                                          .seed = qchem::ChargeDensity::SeedStrategy::SAD}), shared);
    ASSERT_TRUE(calc.Result());
    EXPECT_NEAR(calc.Energy(), -24.40769, 1e-5) << "the shared state, not the integer one (-24.525932)";
}

// D-DEGEN-FILL x GDM (user ruling 2026-11-08, item 3): under DFT the facade's default SHARES boron's partly filled 2p level (one up
// electron, 1/3 per orbital), so D' is not idempotent and GDM must decline it.  The Ladder then stays on its DIIS rung; standalone
// GDM's every step becomes a mixed fixed-point step, and its linear mixer runs under the declined ceiling (0.25, CLAUDE.md: never
// the 1.0 passthrough).  Both reach the SAME shared state -- the energy is the claim that the decline is harmless.
TEST(M_Calculation, BoronLDASharedLevelGDMDeclines)
{
    // The facade converges in its constructor, so the capture brackets the construction.
    testing::internal::CaptureStdout();
    MoleculeCalculation ladder(MoleculeSetup::Resolve(MakeBoron(), {.basis = "dzvp", .model = qchem::Model::LDA, .multiplicity = 2}, {.type = "Ladder"}));
    const std::string outL=testing::internal::GetCapturedStdout();
    testing::internal::CaptureStdout();
    // GDM with the fill UNSTATED takes the integer fill (user 2026-11-08: guide the user); the shared state is stated here.
    SCFParams shared; shared.DegenerateLevels=SCFParams::DegenerateLevelFill::Shared;
    MoleculeCalculation gdm   (MoleculeSetup::Resolve(MakeBoron(), {.basis = "dzvp", .model = qchem::Model::LDA, .multiplicity = 2}, {.type = "GDM"}),
                               shared);
    const std::string outG=testing::internal::GetCapturedStdout();
    ASSERT_TRUE(ladder.Result()) << outL;
    ASSERT_TRUE(gdm.Result())    << outG;
    EXPECT_NE(outL.find("is NOT ENGAGEABLE"), std::string::npos) << "the Ladder must veto GDM on a shared (fractional) D'";
    EXPECT_NE(outG.find("linear mixing capped at relax 0.25"), std::string::npos) << "a declined GDM's fallback runs under the ceiling";
    EXPECT_NEAR(ladder.Energy(), gdm.Energy(), 1e-5) << "one shared state, whichever accelerator declined (to the run's eps)";
}

// Unpolarized (folded-doublet D'=2P) GDM: closed-shell water under standalone GDM must reach the DIIS
// anchor.  GDM used to decline every unpolarized run and silently diagonalize (same energy, no GDM) -- so
// the energy alone cannot tell; the console must carry no "DECLINING" line (checked by hand via ITMain -v).
TEST(M_Calculation, WaterPureGDM)
{
    MoleculeCalculation calc(MoleculeSetup::Resolve(qchem::StructureData::GetMolecule("H2O"), {.basis = "dzvp"}, {.type = "GDM"}));

    const double E_ref = -76.022903;                       // same anchor as WaterSymmetry
    EXPECT_LT(std::fabs((E_ref - calc.Energy()) / E_ref), 1e-5);
}

//====================================================================================================
//  SAVE / RESTART (V-FACADE-STATE, lane A 2026-11-09)
//====================================================================================================
static std::string StatePath(const std::string& leaf)
{
    return (std::filesystem::path(testing::TempDir())/(std::to_string(::getpid())+"_"+leaf)).string();
}

//! A restart's iterates, as its observer sees them (the first one is the claim: CLAUDE.md, pin a restart's FIRST iterate).
struct Iterates
{
    std::vector<double> E;
    MoleculeCalculation::Observer Sink() {return [this](const qchem::SCFIterator::IterationRow& r){ E.push_back(r.eb.GetTotalEnergy()); };}
};

// CLAIM: an EXACT resume starts AT the saved answer -- the saved D drives the first Fock (SeedRole::ExactResume) and the Ladder
// starts on its GDM rung -- so its first iterate is the saved energy, as the NaF gates show for the solid.  Boron UHF/dzvp is the
// open shell (a degenerate 2p level filled by the integer rule).
TEST(M_Calculation, BoronUHFExactResumeStartsAtTheSavedAnswer)
{
    const std::string path=StatePath("boron_uhf.h5");
    const auto setup=[]{ return MoleculeSetup::Resolve(MakeBoron(), {.basis = "dzvp", .multiplicity = 2}, {.type = "Ladder"}); };
    MoleculeCalculation calc(setup());
    ASSERT_TRUE(calc.Result());
    calc.SaveState(path);
    Iterates it;
    auto r=MoleculeCalculation::Restart(path, setup(), {}, it.Sink());
    std::filesystem::remove(path);
    ASSERT_TRUE(r) << r.Error().details;
    ASSERT_TRUE((*r)->RestartedFrom());
    EXPECT_TRUE((*r)->RestartedFrom()->exact);
    ASSERT_TRUE((*r)->Result());
    ASSERT_FALSE(it.E.empty());
    EXPECT_NEAR(it.E.front(), calc.Energy(), 1e-9) << "the first iterate left the saved answer";
    EXPECT_NEAR((*r)->Energy(), calc.Energy(), 1e-9);
    EXPECT_LE((*r)->IterationCount(), 2u);
}

// The same claim for a closed shell under DFT on a symmetry-blocked (SALC) basis: H2O LDA/dzvp, C2v.
TEST(M_Calculation, WaterLDASymmetryExactResumeStartsAtTheSavedAnswer)
{
    const std::string path=StatePath("h2o_lda_c2v.h5");
    const auto setup=[]{ return MoleculeSetup::Resolve(qchem::StructureData::GetMolecule("H2O"),
                                                       {.basis = "dzvp", .model = qchem::Model::LDA, .symmetry = true}); };
    MoleculeCalculation calc(setup());
    ASSERT_TRUE(calc.Result());
    calc.SaveState(path);
    Iterates it;
    auto r=MoleculeCalculation::Restart(path, setup(), {}, it.Sink());
    std::filesystem::remove(path);
    ASSERT_TRUE(r) << r.Error().details;
    ASSERT_TRUE((*r)->RestartedFrom());
    EXPECT_TRUE((*r)->RestartedFrom()->exact) << "the fresh setup must adopt the saved state's XC mesh";
    ASSERT_TRUE((*r)->Result());
    ASSERT_FALSE(it.E.empty());
    EXPECT_NEAR(it.E.front(), calc.Energy(), 1e-9) << "the first iterate left the saved answer";
    EXPECT_NEAR((*r)->Energy(), calc.Energy(), 1e-9);
    EXPECT_LE((*r)->IterationCount(), 2u);
}

// CLAIM (controller note 1, 2026-11-09): POSITIONS ARE WARM, and the density crosses the overlap change through the Löwdin frame.
// A restart at a displaced H2O holds exactly N electrons on the NEW basis, Tr(D S') = 10 to 1e-12, at its first iterate, and
// converges to the fresh run's energy.  DISCRIMINATES (lane A 2026-11-09, the carry disabled): the verbatim AO carry holds
// 9.99506086599 electrons at this 0.064 bohr step.
TEST(M_Calculation, WaterRestartAtMovedAtomsCarriesTheDensityThroughTheLowdinFrame)
{
    const std::string path=StatePath("h2o_moved.h5");
    const Molecule water=qchem::StructureData::GetMolecule("H2O");
    MoleculeCalculation calc(MoleculeSetup::Resolve(water, {.basis = "dzvp"}));
    ASSERT_TRUE(calc.Result());
    calc.SaveState(path);

    std::vector<Vector3D<double>> R;
    water.ForEachSite([&](int, const Vector3D<double>& r, bool){ R.push_back(r); });
    ASSERT_EQ(R.size(), 3u);
    R[1]=R[1]+Vector3D<double>(0.0, 0.05, -0.04);   // one hydrogen, 0.064 bohr
    std::shared_ptr<qchem::Structure> moved=water.Clone();
    moved->MoveSites(R);
    const auto setupMoved=MoleculeSetup::Resolve(*moved, {.basis = "dzvp"});

    auto r=MoleculeCalculation::PrepareRestart(path, setupMoved);
    std::filesystem::remove(path);
    ASSERT_TRUE(r) << r.Error().details;
    const auto how=(*r)->RestartedFrom();
    ASSERT_TRUE(how);
    EXPECT_FALSE(how->exact);
    EXPECT_TRUE(how->lowdinCarried);
    EXPECT_NEAR(how->maxDisplacement, std::sqrt(0.05*0.05+0.04*0.04), 1e-12);
    EXPECT_NEAR(how->electrons, 10.0, 1e-12) << "Tr(D S') on the moved basis";
    bool named=false;
    for (const auto& w : how->warm) named |= w.rfind("structure.positions:",0)==0;
    EXPECT_TRUE(named) << "the warm list names the move";

    auto out=(*r)->Converge();
    ASSERT_TRUE(out) << out.Error().details;
    MoleculeCalculation fresh(setupMoved);
    ASSERT_TRUE(fresh.Result());
    EXPECT_NEAR((*r)->Energy(), fresh.Energy(), 1e-8) << "the carried start reaches the same ground state";
}

// The refusals and the warm differences a binding reads back (JudgeState), one per class: another BASIS is refused (the saved D is
// a matrix over other functions), another MODEL is a warm start (named), and a molecule's state offered to an atom run is refused
// by KIND.
TEST(M_Calculation, ASavedStateIsJudgedByClass)
{
    const std::string path=StatePath("h2o_judged.h5");
    const Molecule water=qchem::StructureData::GetMolecule("H2O");
    MoleculeCalculation calc(MoleculeSetup::Resolve(water, {.basis = "dzvp"}));
    ASSERT_TRUE(calc.Result());
    calc.SaveState(path);

    auto same=calc.JudgeState(path);
    ASSERT_TRUE(same) << same.Error().details;
    EXPECT_TRUE(same->empty()) << "the run that saved it judges it an exact resume";

    auto otherBasis=MoleculeCalculation::Restart(path, MoleculeSetup::Resolve(water, {.basis = "dzvp2"}));
    ASSERT_FALSE(otherBasis);
    EXPECT_EQ(otherBasis.Error().why, qchem::RestartRefusal::Why::Mismatch);
    EXPECT_NE(otherBasis.Error().details.find("basis.name"), std::string::npos) << otherBasis.Error().details;

    auto otherModel=MoleculeCalculation::PrepareRestart(path, MoleculeSetup::Resolve(water, {.basis = "dzvp", .model = qchem::Model::LDA}));
    ASSERT_TRUE(otherModel) << otherModel.Error().details;
    ASSERT_TRUE((*otherModel)->RestartedFrom());
    bool named=false;
    for (const auto& w : (*otherModel)->RestartedFrom()->warm) named |= w.rfind("functional:",0)==0;
    EXPECT_TRUE(named) << "an HF state seeding an LDA run is a WARM start that names the functional";

    auto atom=qchem::AtomCalculation::Restart(path, 8);
    std::filesystem::remove(path);
    ASSERT_FALSE(atom);
    EXPECT_EQ(atom.Error().why, qchem::RestartRefusal::Why::Kind);
}
