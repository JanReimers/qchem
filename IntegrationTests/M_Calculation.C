// File: UnitTests/M_Calculation.C  Tests for the qchem::Calculation front-door facade.
//
// The facade is meant to BE the production recipe -- so its job here is to reproduce, through the
// public front door alone, the converged HF/dzvp water energy that M_HF_U checks via the QchemTester
// scaffold.  Same geometry, same basis, same anchor.  Also a quick smoke of the ScalarFunction
// sampling surface (Density/HOMO/Orbital), the win the API-ergonomics review highlighted.
#include "gtest/gtest.h"
#include <cmath>

import qchem.Calculation;
import qchem.Structure;
import qchem.Types;        // Vector3D
import qchem.Reporting;    // report:: -- the run report the facade emits (scf section)
using namespace qchem;

using qchem::Calculation;

static Molecule MakeWater()      // experimental geometry in BOHR, C2 axis along z (== M_HF_U)
{
    Molecule w;
    w.Insert(new Atom(8, 0, Vector3D<double>(0,  0.0,   0.0)));
    w.Insert(new Atom(1, 0, Vector3D<double>(0,  1.431, 1.107)));
    w.Insert(new Atom(1, 0, Vector3D<double>(0, -1.431, 1.107)));
    return w;
}

// Regression anchor: identical converged HF/dzvp total energy to M_HF_U_Water (-76.022903).
TEST(M_Calculation, WaterEnergy)
{
    Molecule water = MakeWater();
    Calculation calc(water, {.basis = "dzvp"});            // build + converge in one line

    const double E_ref = -76.022903;
    EXPECT_LT(std::fabs((E_ref - calc.Energy()) / E_ref), 1e-5);
    EXPECT_GT(calc.IterationCount(), 0u);
}

// The scf run-report SCHEMA CHECK (RunReportPlan step 1).  The provider (Calculation::Converge)
// builds the `scf` section json directly; here we drive one Converge and assert the emitted shape --
// this is the typo-catcher that replaces the compile-time key-safety we traded for "the report IS json".
TEST(M_Calculation, ScfReportSchema)
{
    Molecule water = MakeWater();
    Calculation calc(water, {.basis = "dzvp"});          // ctor's Converge emits scf into GlobalReport

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
    Calculation calc(MakeWater(), { .basis = "dzvp", .symmetry = true });

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
    Calculation calc(MakeWater(), { .basis = "dzvp" });  // HF over dzvp populates Jac/Kab + the Omega/Hermite caches

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
    Calculation calc(MakeWater(), {.basis = "dzvp", .model = qchem::Model::LDA});
    // Converged parameter-free LDA total energy (Dirac exchange + VWN5 correlation), with the facade's
    // auto DIIS-from-start + SAD seed.  Bit-stable at -75.9324615507 and now ORDER-INDEPENDENT: the
    // ~585 ppm HF-before-DFT drift was a fit-basis Normalization cached without a mesh key (the HF SAD
    // bootstrap's coarse seed mesh poisoned the production run); fixed by mesh-keying Fit_IBS::Norm().
    // So the anchor can stay tight (1e-5) -- if it moves, something regressed.
    const double E_ref = -75.9324615507;
    EXPECT_LT(std::fabs((E_ref - calc.Energy()) / E_ref), 1e-5);
    EXPECT_GT(calc.IterationCount(), 0u);
}

// {.symmetry=true} SALC-blocks the (Cartesian PG) basis + does global aufbau across irreps.  The
// converged total energy must match the un-blocked run -- the M_Sym invariant, now through the facade.
TEST(M_Calculation, WaterSymmetry)
{
    Calculation calc(MakeWater(), {.basis = "dzvp", .symmetry = true});

    const double E_ref = -76.022903;                       // identical to the un-blocked WaterEnergy anchor
    EXPECT_LT(std::fabs((E_ref - calc.Energy()) / E_ref), 1e-5);
    EXPECT_GT(calc.IterationCount(), 0u);
}

// engine=LibCint + {.symmetry=true}: the libcint-Cartesian orbital IBS is ALSO symmetry-adaptable (it
// implements GetAoShells via Molecule::Orbital_1E_IBS, like the MnD bases), so this must converge to the
// same energy as the MnD/un-blocked run.  Guards the {engine=libcint}x{symmetry} combination, previously
// untested -- a symmetry-adaptation refactor once silently regressed it (libcint adapted via the old PGData
// cast); this test locks it in.
TEST(M_Calculation, WaterSymmetryLibCint)
{
    Calculation calc(MakeWater(), {.basis = "dzvp", .engine = Engine::LibCint, .symmetry = true});

    const double E_ref = -76.022903;                       // same anchor as WaterSymmetry (MnD)
    EXPECT_LT(std::fabs((E_ref - calc.Energy()) / E_ref), 1e-5);
    EXPECT_GT(calc.IterationCount(), 0u);
}

// The caller's Molecule is deep-copied: it survives being passed in and may be used afterwards.
TEST(M_Calculation, OwnsItsOwnStructure)
{
    Calculation calc(MakeWater(), {.basis = "dzvp"});      // temporary destroyed at the semicolon
    EXPECT_EQ(calc.GetStructure().GetNumAtoms(), 3u);
    EXPECT_NEAR(calc.GetStructure().GetNumElectrons(), 10.0, 1e-12);
}

// Density and every occupied MO are sampleable ScalarFunction<double>s through one interface.
TEST(M_Calculation, SamplingSurface)
{
    Calculation calc(MakeWater(), {.basis = "dzvp"});

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
    b.Insert(new Atom(5, 0, Vector3D<double>(0, 0, 0)));
    return b;
}

TEST(M_Calculation, BoronUHFLadderGDM)
{
    Calculation calc(MakeBoron(), {.basis = "dzvp", .multiplicity = 2}, {.type = "Ladder"});

    const double E_ref = -24.525932;   // pinned converged UHF/dzvp did-E-move anchor (measured 2026-08-17;
                                       // ~3 mHa above the HF limit -24.529061 = dzvp basis incompleteness)
    EXPECT_LT(std::fabs((E_ref - calc.Energy()) / std::fabs(E_ref)), 1e-5);
    EXPECT_GT(calc.IterationCount(), 0u);
}

// The same anchor under standalone GDM (no DIIS rung in front): the direct minimizer alone must land
// on the same UHF energy.  Locks the second previously-uncovered accelerator type.
TEST(M_Calculation, BoronUHFPureGDM)
{
    Calculation calc(MakeBoron(), {.basis = "dzvp", .multiplicity = 2}, {.type = "GDM"});

    const double E_ref = -24.525932;   // same anchor as the Ladder run (the two agree to ~1e-12 Ha)
    EXPECT_LT(std::fabs((E_ref - calc.Energy()) / std::fabs(E_ref)), 1e-5);
    EXPECT_GT(calc.IterationCount(), 0u);
}
