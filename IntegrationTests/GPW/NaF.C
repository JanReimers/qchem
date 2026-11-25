// File: IntegrationTests/GPW/NaF.C  NaF rocksalt -- the sharp-F ionic multi-species cell (a=8.73, SR2 basis); the best CP2K-validated material (0.2 mHa).
//
// THE GRID (doc/TestSuitePlan.md §3): TEST(<Basis>_<Material>, <k>_[Grid]_[Fit]_[Sym]_[Spin]_[Occ]_[Machinery]_[Ansatz]_[Seed]_<Claim>)
// -- axis tokens in fixed order, the facade's defaults elided (k always named); the claim is CP2K (oracle anchor),
// Anchor (did-E-move), eq<Token> (a ONE-axis twin: this point equals the same point with that axis moved) or a
// property verb.  Tests are laid out in axis order.  `scripts/testgrid` renders the coverage table from these names.
//
//   GPW_NaF.Γ_Imp_Anchor
//   GPW_NaF.Γ_Becke_Imp_eqUni
//   GPW_NaF.Γ_Imp_eqColdStart       (coarse -> Restart onto the fine grid == a cold fine run)
//   GPW_NaF.Γ_Imp_eqExactResume     (same-grid Restart reproduces the converged state)
//   GPW_NaF.Γ_PathIsContinuous_Long (a SolidSession's 20-point displacement path has no jump: CalcNetworkReview gate (c))

#include "gtest/gtest.h"
#include <unistd.h>    // getpid: per-PROCESS temp names (ctest runs each case as its own process, in parallel)
#include <memory>
#include <vector>
#include <cmath>
#include <cstdlib>   // std::getenv/std::atof (the NaF mixing-tuning env knobs)
#include <complex>
#include <cstdio>
#include <fstream>   // /proc/self/statm (the RSS breadcrumb bisect)
#include <stdexcept>
#include <algorithm>
#include <functional>
#include <string>
#include <filesystem>
#include <iostream>
#include <iomanip>      // setprecision (the order-parameter trajectory line)
#include <nlohmann/json.hpp>   // the eps gate's inline decks

import qchem.SCFCap;                             // IterationBudget: the convergence PANEL's budget (doc/SCFRobustnessPlan.md §4)
import qchem.UnitCell;                           // UnitCell, FCCUnitCell
import qchem.Lattice_3D;                         // Lattice_3D
import qchem.BasisSet;                           // Complex_BS, Real_BS
import qchem.BasisSet.Orbital_1E_IBS;            // Complex_OIBS (the overlap-spectrum diagnostic)
import qchem.BasisSet.PlaneWave.PlaneWave_IBS;   // PlaneWave_IBS (the seed's CD fit basis)
import qchem.BasisSet.Lattice.BasisSet;       // GPWFactory (the GPW basis container)
import qchem.BasisSet.Gaussian.Point.Factory;          // Gaussian::Factory, BasisSetData/Engine/Angular
import qchem.Outcome;                           // Outcome<Converged,SCFFailure> -- the facade's result
import qchem.SolidCalculation;                    // the NAMED periodic facade (Step 4 3/3)
import qchem.Tests.GPW_Harness;                   // THE HARNESS (IntegrationTests/GPW/Harness.C): Materials cells, gates, recipes, the XC probes
import qchem.Materials;                           // Materials::Get -- the cells come from src/Structure/Data/materials.json (row MD)
import qchem.ChargeDensity.DensitySampler;  // the XC sampling engine (its own module
                                                  // since 2026-09-08; .Internal. modules are
                                                  // never re-exported, so name it directly)
import qchem.BasisSet.G_FieldEvaluator;           // G_RasterTransform -- the uniform probe's own point count
import qchem.Mesh.Angular;                        // MakeAngular (the rotated-Lebedev bond-angle probe)
import qchem.Mesh;                                // qcMesh::MeshParams / UnitCellKind (the Becke XC quadrature)
import qchem.Mesh.XCPolicy;                       // BeckeXCParams / ResolveXCMesh / XCMeshSharpness (the grid policy)
import qchem.SCFIterator;                        // cSCFIterator, SCFParams
import qchem.SCFParams;                          // SCFParams
import qchem.ChargeDensity.Seed;                 // SeedStrategy
import qchem.SCFAccelerator.Factory;              // the PUBLIC complex accelerator door (Step 4)
import qchem.WaveFunction;                       // cWaveFunction (the converged state)
import qchem.Energy;                             // EnergyBreakdown
import qchem.Symmetry.Irrep;                     // Irrep
import qchem.Symmetry.Spin;                      // Spin
import qchem.LASolver;                           // qchem::Ortho (Cholesky | Eigen | SVD -- basis orthogonalisation)
import qchem.BasisSet.Gaussian.Lattice.GPW_IBS;         // GPW_IBS (build a concrete block for the collocation diagnostic)
import qchem.BasisSet.GMap;              // PairProjection<dcmplx> (the collocation weight tensor); SymmetryDefects (§3 diagnostic)
import qchem.Deck;                                // deck::Run -- the eps gate runs the path rundeck and the GUI run
import qchem.Reporting;                           // the eps gate reads the run's XC-mesh policy row
import qchem.Pseudopotential.GTH_Potentials;      // GetGTH, GTH_PP (the PP model, for the matrix-trace probe)
import qchem.Types;
import qchem.GeometryOptimizer;                   // GeometryOptimizer, CellWyckoffMap (the rocksalt cell relaxation)
import qchem.Matrix3D;                           // Matrix3D<double> (the fcc start cell)

using namespace qchem;
using BasisSet::Real_BS;
using BasisSet::Complex_BS;
using qchem::BasisSet::Gaussian::BasisSetData;

using namespace qchem::tests::gpw;   // the harness: Materials cells, gates, recipes, probes (IntegrationTests/GPW/Harness.C)


// (4) MULTI-SPECIES GPW: ionic NaF (rocksalt = FCC + 2-atom basis) at Gamma, driven by the multi-species
// Ham_PW_DFT ctor ({{"Na",1},{"F",7}}).
// Valance basis for Na generated from all electrton atom caluclations (BasisSetData/valence_lowq.bsd) 
// looks like     s={0.03.0.086,0.245,0.7,2.0}, p={0.05, 0.3}
// The diffuse exponents s={0.03.0.086} simply don't work in a lattice context.  The overlap is singular.  Auto trimming
// the basis using eeigen instead choleslky decomposition also simply doe not work.  The only option is to tell the
// user to drop thos diffuse basis function.  If use the BasisSetData::VALENCE_LOWQ_SR2 the calculation converges nicely
// with some help from DIIR/GDM/Kerker-mixing.
// Also work noting is the sharp F basis function with exponent 40Ha. THis motivated a ladder of grids to which basis 
// function pairs are assigned bases on thier exponents.  
// There is also a challange for the integration grids used for Vxc fitting.  It is anticipated that using a non-uniform unit Becke/Voronoi-polyhedra grid will allow
// for rapid integration of diffuse basis function will make this run even more efficient. 
// The ideal minimum densityEcut=2*40Ha=80Ha based on the F max exponent.  40 converges to a lower E_total but otherwise converged nicely.
// Explicit densityEcut= 40 = SUB-FLOOR: warns, and the ball's own raster aliases there (-43 mHa)
TEST(GPW_NaF, Γ_Imp_Anchor)   // RE-ENABLED 2026-09-15: 15 s, converged in 23 iterations -- the GPW x NaF row's first standing anchor
{
    // THE NaF ORACLE ANCHOR (doc/GPWPlan.md: 0.10-0.19 mHa vs CP2K on this basis, tight-eps + converged density).
    // The committed production recipe (NaFOptions/NaFGates), SR2 basis, Γ -- the arm doc/Benchmark.md times and
    // the CP2K decks (naf_gpw_sr2_diag.inp / naf_gpw_sr_tight.inp, no &KPOINTS) compute.  The k-mesh, span,
    // ecut, ladder, alpha, smearing and penalty knobs this test used to read from NAF_* env vars are the
    // campaign INSTRUMENT's and go with it to the probe binary (doc/TestSuitePlan.md §8).
    const Material naf=qchem::Materials::Get("NaF_rocksalt");
    const Lattice_3D lat=LatticeOf(naf);
    RunOptions o=NaFOptions(naf, "NaF GPW Gamma");
    o.imposeSymmetry=true;   // V1.30: was the DEFAULT; now stated, because an imposition you did not ask for is invisible in the result
    o.xcMesh          = qcMesh::BeckeXCParams(20,2,24);   // BELOW the calibrated default (nR 40, degree 29): this anchor's number only, never a curvature (OpenWork "A curvature needs an XC-MESH ladder")
    SolidXCMesh(o).cellKind=qcMesh::CellKind::RecipeStated(qcMesh::UnitCellKind::Becke);
    SCFParams par=NaFGates(); par.StartingRelaxRo=0.25; par.MergeTol=1e-4;
    GpwReport report("NaF "+o.label, par.Verbose);
    o.scf=par; auto calc=qchem::SolidCalculation::Prepare(SolidSetup::Resolve(lat, MakeBasisNaFSR2(*naf.cell), o)).TakeValue(); (void)calc->Converge();
    auto R=calc->Result();
    ASSERT_TRUE(R) << Why(R);
    EXPECT_NEAR(R->TotalCharge(), 8.0, 1e-6);           // 1 (Na) + 7 (F) valence electrons, conserved
    EXPECT_NEAR(R->Energy(), -24.4304, 0.01);           // did-E-move anchor at Γ (CP2K agrees to 0.2 mHa, doc/GPWPlan.md)
}


// INCREMENTAL CONVERGENCE ON A HARD MATERIAL (D-NAFGRID, user 2026-10-02: "very important ... well tested, and on a
// hard material like NaF").  NaF is the stress case for a restart: the sharp F pseudopotential and the diffuse Na
// valence basis make it easy to land in the wrong state after a perturbation -- the old raw-iterator version of this
// test (GPW_NaF.DISABLED_Γ_GridContinuation) recorded a transferred seed PINNING AN EXCITED STATE (-23.68 against
// -24.43 Ha).  So: converge a CHEAP coarse density grid, save, and Restart onto the production grid.  Same orbital
// basis, same recipe as Γ_Imp_Anchor; only densityEcut (and, for the coarse stage, the alias-free raster the sub-floor
// Ecut=40 needs) differs.  One fixture, built once: the coarse state, the cold fine reference, and the fine state.
namespace
{
struct NaFRestartFixture
{
    std::string coarsePath, finePath;
    double      coldEnergy=0.0, coarseEnergy=0.0;
    size_t      coldIterations=0;
    std::shared_ptr<const SolidSetup> fineSetup;   //!< the cold fine run's setup, its rasters as its pilot sized them
    bool        ok=false;
    std::string why;
};
RunOptions NaFAnchorOptions(const Material& naf, const std::string& label)
{
    RunOptions o=NaFOptions(naf, label);
    o.imposeSymmetry=true;
    o.xcMesh          = qcMesh::BeckeXCParams(20,2,24);   // BELOW the calibrated default (nR 40, degree 29): this anchor's number only, never a curvature (OpenWork "A curvature needs an XC-MESH ladder")
    SolidXCMesh(o).cellKind=qcMesh::CellKind::RecipeStated(qcMesh::UnitCellKind::Becke);
    return o;
}
SCFParams NaFAnchorParams() { SCFParams par=NaFGates(); par.StartingRelaxRo=0.25; par.MergeTol=1e-4; return par; }

const NaFRestartFixture& NaFRestart()
{
    static const NaFRestartFixture f=[]
    {
        NaFRestartFixture r;
        const Material naf=qchem::Materials::Get("NaF_rocksalt");
        const Lattice_3D lat=LatticeOf(naf);
        const auto dir=std::filesystem::path(testing::TempDir());
        // PER PROCESS: ctest runs Γ_Imp_eqColdStart and Γ_Imp_eqExactResume as separate processes, in parallel, each
        // building this fixture -- a shared name let one restart from the other's half-written file (2026-10-21).
        const std::string pid=std::to_string(::getpid());
        r.coarsePath=(dir/("naf_restart_coarse."+pid+".h5")).string();
        r.finePath  =(dir/("naf_restart_fine."+pid+".h5")).string();
        {   // the CHEAP coarse stage: Ecut=40 is SUB-FLOOR (below C*alpha_max=80), so it needs a finer raster than its ball's own --
            // stated (accuracy 0): x 16, twice the ball's G_N (C2 step 3)
            RunOptions o=NaFAnchorOptions(naf, "NaF restart coarse");
            o.solid.densityEcut=40.0;
            o.accuracy.eps=0.0;
            o.solid.tolerances.rasterX=o.solid.tolerances.vlocRasterX=o.solid.tolerances.xcRasterX=16.0;
            o.saveStateTo=r.coarsePath;
            o.scf=NaFAnchorParams(); auto calc=qchem::SolidCalculation::Prepare(SolidSetup::Resolve(lat, MakeBasisNaFSR2(*naf.cell), o)).TakeValue(); (void)calc->Converge();
            auto R=calc->Result();
            if (!R) { r.why="coarse: "+Why(R); return r; }
            r.coarseEnergy=R->Energy();
        }
        {   // the COLD fine reference (what Γ_Imp_Anchor runs), also saved for the exact-resume claim
            RunOptions o=NaFAnchorOptions(naf, "NaF restart cold fine");
            o.saveStateTo=r.finePath;
            o.scf=NaFAnchorParams(); auto calc=qchem::SolidCalculation::Prepare(SolidSetup::Resolve(lat, MakeBasisNaFSR2(*naf.cell), o)).TakeValue(); (void)calc->Converge();
            auto R=calc->Result();
            if (!R) { r.why="cold fine: "+Why(R); return r; }
            r.coldEnergy=R->Energy(); r.coldIterations=R->IterationCount();
            r.fineSetup=calc->SizedSetup();   // the PRODUCTION grid, as this run's pilot sized it (C2 step 3): what the coarse state restarts onto
        }
        r.ok=true;
        return r;
    }();
    return f;
}
} // namespace

// CLAIM: converge coarse, Restart onto the fine grid, and land on the SAME ground state a cold fine run reaches --
// not an excited-state pin, not the -40 Ha basin.  The first iterate is judged too: the whole point of the seed is
// that the fine SCF STARTS near the answer (a count of iterations is the accelerator's business, SolidState.C).
TEST(GPW_NaF, Γ_Imp_eqColdStart_Panel)   // PANEL case 2 (regime 5, WARM coarse -> fine): the fine stage's budget
{
    const NaFRestartFixture& f=NaFRestart();
    ASSERT_TRUE(f.ok) << f.why;
    qchem::IterationBudget budget(5, "NaF coarse -> fine restart, the fine stage (regime 5; 3 and 4 measured 2026-10-29)");
    EXPECT_NEAR(f.coarseEnergy, -24.4357, 0.01) << "the coarse seed stage must itself be the physical fixed point";
    const Material naf=qchem::Materials::Get("NaF_rocksalt");
    const Lattice_3D lat=LatticeOf(naf);
    RunOptions o=NaFAnchorOptions(naf, "NaF restart fine");          // auto densityEcut = the production grid
    ASSERT_EQ(o.accelerator, qchem::SCFAccelerators::Type::Ladder);        // the production recipe; Restart enters at GDM
    std::vector<double> E;
    // Onto the cold fine run's SIZED setup (C2 step 3): its rasters are the production grid's, so the claim -- the first fine iterate is
    // already the answer -- is about the restart, not about a pilot re-sizing from the ball's rasters.
    ASSERT_TRUE(f.fineSetup);
    auto c=qchem::SolidCalculation::Restart(f.coarsePath, f.fineSetup, NaFAnchorParams(),
                                            [&E](const qchem::SCFIterator::IterationRow& p){ E.push_back(p.eb.GetTotalEnergy()); });
    ASSERT_TRUE(c) << (c ? std::string() : c.Error().details);
    auto R=(*c)->Result();
    ASSERT_TRUE(R) << Why(R);
    ASSERT_FALSE(E.empty());
    EXPECT_NEAR(R->TotalCharge(), 8.0, 1e-6);
    EXPECT_NEAR(R->Energy(), f.coldEnergy, 1e-4) << "the restarted run reached a different state than the cold fine run";
    EXPECT_NEAR(R->Energy(), -24.4304, 0.01) << "and it is the banked NaF anchor";
    // THE LADDER ENTERS AT ITS GDM RUNG on a restart (user 2026-10-21).  MEASURED 2026-10-02 before that: from the
    // coarse density the fine run STARTS 4e-8 Ha from the answer and GDM converges in 2 iterations, while the Ladder's
    // DIIS+Kerker rung stepped away (1.5e-4) and took 12 (cold fine run: 21).  A restart whose density moved MORE
    // (geometry, U) is unmeasured.
    // RE-MEASURED 2026-10-23 after the V_short operator-image fix (BuildOpTerms): the first fine iterate is now 2.3e-6 Ha from the
    // cold fine answer (was 4e-8) -- NOT convergence slop (seeds converged to Δρ/N 1e-6 gave the same 2.27e-6); the coarse grid's
    // fixed point sits further from the fine one than before.  Still 65x closer than the DIIS rung's 1.5e-4, and the final energy
    // matches the cold run to 1e-8.  Why the gap grew is open (OpenWork §3 "MnO accuracy", the V_short fix entry).
    EXPECT_LT(std::abs(E.front()-f.coldEnergy), 1e-5) << "the first fine iterate should already be (nearly) the answer";
    EXPECT_LE(R->IterationCount(), 4u) << "the coarse seed must save iterations on the fine grid (measured 2; cold fine run 21)";
}

// CLAIM: the SAME grid resumes EXACTLY -- Restart from the converged fine state reproduces its energy and needs
// (almost) no iterations.  The control for the claim above: if this fails the file round-trip is lossy on a hard
// material, and the grid-continuation result means nothing.
TEST(GPW_NaF, Γ_Imp_eqExactResume_Panel)   // PANEL case 1 (regime 4, exact resume)
{
    const NaFRestartFixture& f=NaFRestart();
    ASSERT_TRUE(f.ok) << f.why;
    qchem::IterationBudget budget(3, "NaF exact resume (regime 4)");
    const Material naf=qchem::Materials::Get("NaF_rocksalt");
    const Lattice_3D lat=LatticeOf(naf);
    RunOptions o=NaFAnchorOptions(naf, "NaF restart exact");
    ASSERT_EQ(o.accelerator, qchem::SCFAccelerators::Type::Ladder);        // the production recipe; Restart enters at GDM
    std::vector<double> E;
    // Onto the cold fine run's SIZED setup, as above: a fresh Resolve is the ball's auto grid, which the saved state is not a fixed point of.
    ASSERT_TRUE(f.fineSetup);
    auto c=qchem::SolidCalculation::Restart(f.finePath, f.fineSetup, NaFAnchorParams(),
                                            [&E](const qchem::SCFIterator::IterationRow& p){ E.push_back(p.eb.GetTotalEnergy()); });
    ASSERT_TRUE(c) << (c ? std::string() : c.Error().details);
    auto R=(*c)->Result();
    ASSERT_TRUE(R) << Why(R);
    ASSERT_FALSE(E.empty());
    // THE LADDER'S GDM RUNG (user 2026-10-21: a restart enters there).  MEASURED 2026-10-02, before that: from the same
    // file the Ladder's DIIS rung + Kerker(alpha=.25) STEPPED the density before the first energy was read -- first
    // iterate 1.3e-4 Ha off the saved state, 4 DIIS iterations until |dE/E|<1e-6 handed off to GDM, 9 in all -- whereas
    // GDM starts at the answer (4e-10) and converges in 2.  The saved state was never the problem.
    EXPECT_NEAR(E.front(), f.coldEnergy, 1e-7) << "the first iterate of an exact resume is the saved state";
    EXPECT_LE(R->IterationCount(), 3u) << "an exact resume on the Ladder's GDM rung starts converged";
    EXPECT_NEAR(R->Energy(), f.coldEnergy, 1e-7);
}


// The SHARP-FIELD leg of the gate (the plan names DISABLED_NaFRocksaltGamma as the stress case: the F-
// anion makes sharp peaks in rho and V_xc, and its diffuse basis is what the Becke grid exists for).
// DISABLED like the parent NaF anchor -- it is a long run (the full NaF convergence recipe at a
// matrix-grade densityEcut=160 reference, plus two Becke term evaluations); run it by hand with
// --gtest_also_run_disabled_tests when touching the XC quadrature.
TEST(GPW_NaF, Γ_Becke_Imp_eqUni)   // RE-ENABLED 2026-09-15: 24 s, Becke internally converged on the sharp-F system
{
    const Material naf=qchem::Materials::Get("NaF_rocksalt");
    const Lattice_3D lat=LatticeOf(naf);

    // The committed NaF production recipe (DISABLED_NaFRocksaltGamma) at PRODUCTION grids (auto Ecut=80,
    // the ball's own raster): the SCF only supplies the density; the comparison itself carries the reference-grade work.
    // MEASURED ATTRIBUTION (2026-07-30): vs the SAME Becke B40 matrix, the uniform reference gave
    //   the ball's own raster, Ecut=160: max|U-B|=1.55e-1      (the production raster's sharp-F element error)
    //   a raster at twice G_N, Ecut=160: max|U-B|=1.89e-2      (8x better -- most of the gap was the ball raster's)
    // with dExc tiny throughout (1.6e-4 / 1.3e-5) -- the discrepant elements carry no rho weight.  So on a
    // sharp-F system the UNIFORM raster is not element-converged at practical Ecut, which is this grid's
    // reason to exist; the gate here is Becke INTERNAL convergence (B40 vs a 2x-refined B80) plus the
    // energy-level agreement with the uniform route.
    RunOptions o=NaFOptions(naf, "NaF Becke gate");
    o.imposeSymmetry=true;   // V1.30: was the DEFAULT; now stated, because an imposition you did not ask for is invisible in the result
    GpwReport report("NaF "+o.label, false);
    o.scf=NaFGates(); auto calc=qchem::SolidCalculation::Prepare(SolidSetup::Resolve(lat, MakeBasisNaFSR2(*naf.cell), o)).TakeValue(); (void)calc->Converge();
    ASSERT_TRUE(calc->Result()) << Why(calc->Result());
    EXPECT_NEAR(calc->Result()->TotalCharge(), 8.0, 1e-6);
    const GpwHandles h=Handles(*calc);

    auto st=lat.GetStructure();
    XCProbe U  =UniformXCProbe(h, st);
    XCProbe B40=BeckeXCProbe(h, st, "B40", qcMesh::BeckeXCParams());
    XCProbe B80=BeckeXCProbe(h, st, "B80", qcMesh::BeckeXCParams(/*nRadial*/80, /*mhlAlpha*/1.0, /*L*/29));
    EXPECT_NEAR(B40.Exc, U.Exc, 2e-3);               // energy-level agreement with the uniform route
    EXPECT_NEAR(B40.rhoLost, 0.0, 5e-3);
    DiffXC(U,B40);                                    // report-only: the uniform raster's element error
    EXPECT_LT(DiffXC(B40,B80), 2e-3) << "Becke V_xc not internally converged on the sharp-F system";
}


// CLAIM (D-BECKE-ENERGY-AUDIT, controller 2026-11-24; user: "I was playing eps and becke grid would not budge"): an LDA ENERGY run on the Becke
// mesh its facade chose (xcMesh unstated) sizes that mesh FROM EPS, and its policy row says so.  The solid now takes the molecule's rule
// (qcMesh::BeckePilotRoleOf: unstated climbs, stated is audited) on every run, not only a geometry run or a GGA; until then an LDA energy run
// kept (40, L29) whatever eps said -- on MnO FM that mesh was 161 uHa low at its own D (OpenWork §3 "Becke vs uniform XC").
// MEASURED 2026-11-24 (rundeck, `~/Code/qchem6-runs/NaF/laneB_becke_eps/`): eps 1e-4 and 1e-5 settle on (40, L29) (radial |dE_xc| 1.44e-5
// against 2e-5 at 1e-5), 5e-6 climbs to (60, L29) in one raise, 1e-6 to (60, GL41) in three.  Two arms here, ~2 + ~3.5 min.
namespace
{
//! The XC-mesh policy rows (SolidXCMeshSettled / NotSettled) of the LAST run filed under \a prefix.
std::vector<qchem::report::json> XCMeshRows(const std::string& prefix)
{
    const auto& G=qchem::report::GlobalReport();
    const qchem::report::json* last=nullptr;
    for (auto it=G.begin(); it!=G.end(); ++it) if (it.key().rfind(prefix,0)==0) last=&it.value();
    if (!last) throw std::logic_error("no filed run "+prefix);
    std::vector<qchem::report::json> rows;
    std::function<void(const qchem::report::json&)> walk=[&](const qchem::report::json& x)
    {
        if (x.is_object())
            for (auto it=x.begin(); it!=x.end(); ++it)
            {
                if (it.key()=="policy" && it.value().is_array())
                    for (const auto& r : it.value()) if (r.contains("name") && r["name"].get<std::string>().rfind("SolidXCMesh",0)==0) rows.push_back(r);
                walk(it.value());
            }
        else if (x.is_array()) for (const auto& y : x) walk(y);
    };
    walk(*last);
    return rows;
}
} // namespace
TEST(GPW_NaF, Γ_Imp_EnergyRunSizesItsBeckeMeshFromEps_Long)
{
    namespace fs=std::filesystem;
    auto run=[](double eps)
    {
        const nlohmann::json deck={
            {"structure","NaF_rocksalt"}, {"kmesh",{1,1,1}}, {"basis",{{"data","VALENCE_LOWQ_SR2"}}}, {"seed","IonicSAD"},
            {"ortho","CholeskyPivoted"}, {"orthoTol",1e-4}, {"accelerator","Ladder"}, {"imposeSymmetry",true}, {"accuracy",eps},   // NO xc mesh: the facade's choice
            {"scf",{{"NMaxIter",100},{"minDeltaRhoPerElectron",1e-7},{"minFD",1e-7},{"minDeltaE","off"},{"minDeltaFD","off"},{"minVirial","off"}}}};
        qchem::deck::RunSpec spec; qchem::deck::FromJson(deck, spec);
        const fs::path dir=fs::temp_directory_path()/("qchem_it_naf_becke_eps_"+std::to_string(::getpid()));
        fs::remove_all(dir);
        qchem::deck::Provenance pv; pv.codeVersion="it";
        const auto out=qchem::deck::Run(spec, pv, dir);
        fs::remove_all(dir);
        EXPECT_TRUE(out.converged) << "eps " << eps << ": " << out.summary;
        const auto rows=XCMeshRows("NaF_rocksalt");
        EXPECT_EQ(rows.size(), 1u) << "eps " << eps << ": one XC-mesh policy row";
        if (rows.empty()) return std::string();
        EXPECT_EQ(rows[0]["name"], "SolidXCMeshSettled") << rows[0].dump();
        EXPECT_EQ(rows[0]["verdict"], "choose");
        const std::string decided=rows[0]["decided"];
        std::cout << "[eps gate] eps " << eps << ": " << decided << std::endl;
        return decided.substr(0, decided.find(" (entered at"));   // the mesh it settled on
    };
    const std::string loose=run(1e-5), tight=run(5e-6);
    EXPECT_NE(loose.find("nR=40, L=29"), std::string::npos) << "eps 1e-5 keeps the entry rung (radial 1.44e-5 against 2e-5): " << loose;
    EXPECT_NE(loose, tight) << "a tighter eps moves the Becke mesh on an LDA energy run";
    EXPECT_EQ(tight.find("nR=40"), std::string::npos) << "eps 5e-6 raises the radial axis (measured (60, L29)): " << tight;
}


// CLAIM (user ruling 2026-11-01, pin 31): a geometry point with NO XC-mesh choice -- xcMesh Auto -- takes the atom-centred Becke mesh and
// its pilot sizes that mesh on the nested ladder (qcMesh::BeckeAngularRungAbove, BeckeRadialRungAbove) to max|dH_xc| <= eps per element; gated by an INDEPENDENT
// audit two rungs above the settled mesh.  The egg-box cannot gate this (the mesh moves with the atoms); the nested ladder is the gate.
// Why the default recipe is not enough here: lane B's NaF/LiF curvature (2026-11-01) -- (20, 2, 24), the path test's mesh below, was 4-10x
// wrong on a 0.1 bohr curvature, and (40, 2, 36) brought both crystals to CP2K per term (NaF dExc -3.08e-4 vs -3.09e-4).  MEASURED 2026-11-01:
// rung differences 1.0e-4 (40/29 vs 60/35 -- the default is 10x over), 1.8e-5 (60/35 vs 80/41 GL), 2.3e-6 (80/41 vs 100/47); settles on
// 80/41 in 3 rounds, 5 re-convergence iterations; against 120/53 it reads 2.7e-6.  ~520 s serial.
TEST(GPW_NaF, Γ_GeometryRunSizesItsBeckeMesh_Long)
{
    const Material naf=qchem::Materials::Get("NaF_rocksalt");
    RunOptions o=NaFOptions(naf, "NaF geometry XC mesh");
    o.imposeSymmetry=false;
    o.solid.xcMeshSymmetric=false;
    o.xcMesh=qcMesh::MeshParams{.cellKind=qcMesh::UnitCellKind::Auto};
    SCFParams par=NaFAnchorParams(); par.MinΔρPerElectron=1e-6;
    const GaussianBasisRecipe sr2{.data=BasisSetData::VALENCE_LOWQ_SR2};
    qchem::SolidSession pes(SolidSetup::Resolve(LatticeOf(naf), sr2, o), {SCFStage{par, o.accelerator}});
    const PESGeometry ref=pes.Reference();
    auto P=pes.ConvergeAt(ref.cell, ref.R);
    ASSERT_TRUE(P) << P.Error().details;
    const SolidCalculation& c=*pes.Current();
    const qcMesh::MeshParams xc=c.ResolvedXCMesh();
    ASSERT_EQ(xc.cellKind, qcMesh::UnitCellKind::Becke) << "a geometry run's XC follows the atoms";
    EXPECT_TRUE(c.PilotCost().settled);
    const qcMesh::BeckeLadderCalibration cal;   // two rungs up on the angular axis, one on the radial
    const auto A=c.XCMeshAudit(*qcMesh::BeckeRadialRungAbove(*qcMesh::BeckeAngularRungAbove(*qcMesh::BeckeAngularRungAbove(xc, cal, false), cal, false)));
    ASSERT_TRUE(A) << A.Error().Summary();
    double dXC=0.0;
    for (const auto& t : A->terms) if (t.entries=="Exc") dXC=std::max(dXC, t.maxAbsΔH);
    std::cout<<"[XC mesh] settled nR="<<xc.nRadial<<" L="<<xc.angularDegree<<"  max|dH_xc| two rungs up "<<dXC<<std::endl;
    EXPECT_GT(dXC, 0.0) << "the audit compared two different meshes";
    EXPECT_LE(dXC, o.accuracy.eps) << "max|dH_xc| against a mesh two rungs finer";
}

// (The pre-facade grid-continuation test, DISABLED_Γ_GridContinuation -- raw SolidSCFIterator ctors, ~15 GC_* env
// knobs -- was DELETED 2026-10-02 (D-NAFGRID): its basin-avoidance premise (the direct fine run diving to -40 Ha) no
// longer reproduces, and its incremental-convergence claim lives in Γ_Imp_eqColdStart / Γ_Imp_eqExactResume above.
// History: doc/Records/OpenWork_History*.md, doc/OldPlans/GPWPlan.md §0e.)

// CLAIM (doc/CalcNetworkReview.md §5 step 8, gate (c)): the potential energy surface is CONTINUOUS along a displacement path.  A
// SolidSession walks F through 20 points of a Γ optical-phonon path, each point seeded from the previous one's density, every point
// on the basis the setup FROZE at the reference (the recipe's trim, never re-decided per step).  A re-trim, a re-folded group, a
// carried density that lands in another basin, or a stale geometry term would each show as a JUMP.  The jump measure is the FIFTH
// difference: on a smooth anharmonic curve the k-th differences fall geometrically with k (measured here: Δ³ 4.5e-7, Δ⁴ 1.0e-7,
// Δ⁵ 2.0e-8 Ha), while a step s at one point shows in Δ⁵ as up to 10 s -- so 1e-7 catches any jump of ~1e-8 Ha.  (A third
// difference conflates the jump with the curve's own E'''h³, which is 4.5e-7 at the path's ends.)  The path is symmetric about the
// site, and rocksalt's mirror through F makes E(+δ) = E(−δ) -- a second invariant that a path walked in ONE direction (carry from
// −δ toward +δ) has to meet.  CAUGHT (2026-10-30, the claim's first run): a free run's XC mesh was built invariant under the group
// DETECTED at each geometry -- 48 ops on the site, 8 off it -- so the mesh switched along the path: mirror 1.7e-5 Ha, Δ⁵ 5.9e-6;
// a SolidSession now requires the plain atom-centred mesh (RunOptions::xcMeshSymmetric=false): 5e-10 and 2.0e-8.
// NOT claimed: that the site is a minimum.  With GTH-q1 Na it is not, in CP2K either (a weak maximum, E'' ≈ -0.005 Ha/bohr² at k222
// on a resolved grid; q9 Na gives +0.029 ours / +0.028 CP2K on the same basis) -- doc/OpenWork.md "A light-PP ion has NO restoring
// force on its site".  (The -0.03 once quoted here was this run's raster egg-box: its rigid-translation control moves E 1e-5.)
TEST(GPW_NaF, Γ_PathIsContinuous_Long)
{
    const Material naf=qchem::Materials::Get("NaF_rocksalt");
    RunOptions o=NaFOptions(naf, "NaF PES path");
    o.imposeSymmetry=false;   // a displaced F lowers the group: an imposed session throws on it: a displaced-atom driver runs a lower-group setup
    o.solid.xcMeshSymmetric=false;  // ... and a free session's XC mesh must not follow the group detected at each point (see above)
    o.xcMesh          = qcMesh::BeckeXCParams(20,2,24);   // BELOW the calibrated default: the CONTINUITY claim only -- a site curvature on it is 4-10x off (LiF/NaF, OpenWork "A curvature needs an XC-MESH ladder")
    SolidXCMesh(o).cellKind=qcMesh::CellKind::RecipeStated(qcMesh::UnitCellKind::Becke);
    o.solid.densityEcut     = 80.0;   // STATED, the old auto ball (2 alpha_max): the pilot would double it (Hartree 1.9e-2 over at 80 Ha) and the 20 serial points cost 2574 s, not ~600 -- the CONTINUITY claim does not need the ball converged, only fixed along the path
    SCFParams par=NaFAnchorParams(); par.MinΔρPerElectron=1e-6;
    const GaussianBasisRecipe sr2{.data=BasisSetData::VALENCE_LOWQ_SR2};
    qchem::SolidSession pes(SolidSetup::Resolve(LatticeOf(naf), sr2, o), {SCFStage{par, o.accelerator}});
    const PESGeometry ref=pes.Reference();
    const auto f=std::find(ref.Z.begin(), ref.Z.end(), 9);
    ASSERT_NE(f, ref.Z.end());
    const size_t iF=f-ref.Z.begin();

    constexpr size_t N=20;
    constexpr double h=0.01;   // bohr, Cartesian x; the path runs δ = -0.095 ... +0.095
    std::vector<double> E(N), d(N);
    std::vector<size_t> it(N);
    for (size_t i=0;i<N;i++)
    {
        d[i]=(double(i)-0.5*(N-1))*h;
        std::vector<rvec3_t> R=ref.R;
        R[iF].x+=d[i];
        auto P=pes.ConvergeAt(ref.cell, R);
        ASSERT_TRUE(P) << "δ=" << d[i] << ": " << P.Error().details;
        E[i]=P->energy; it[i]=P->iterations;
        std::cout << "[PES path] δ=" << std::setw(7) << d[i] << "  E=" << std::setprecision(12) << E[i] << "  iterations=" << it[i] << std::endl;
    }
    for (size_t i=0;i<N/2;i++)
        EXPECT_NEAR(E[i], E[N-1-i], 1e-7) << "the mirror through F: E(δ) == E(-δ) at δ=" << d[N-1-i];
    for (size_t i=0;i+5<N;i++)
    {
        const double d5=E[i+5]-5*E[i+4]+10*E[i+3]-10*E[i+2]+5*E[i+1]-E[i];
        EXPECT_LT(std::abs(d5), 1e-7) << "a jump between δ=" << d[i] << " and δ=" << d[i+5];
    }
}

//------------------------------------------------------------------------------------------------------------------------------------
// ROCKSALT NaF RELAXES ITS CELL ON GPW (doc/OpenWork.md §2 "Real relaxations end to end", the easy cases: controller plan 2026-11-20 (b)).
// Cell only: rocksalt has no free internal coordinate, so the map over the detected group (Fm-3m) has ONE coordinate, the cubic constant.
// FREE density (all 8 k-points, the plain Becke mesh), the group held by the map.  SEMICORE Na q9: with q1 Na the site has no restoring
// force, in CP2K too (OpenWork "A light-PP ion has NO restoring force"), so a q1 crystal is no stable structure to relax.  The oracle is
// CP2K CELL_OPT on the identical model (`~/Code/cp2k-runs/NaF/cellopt/naf_q9_k222_cellopt.{inp,out}`): the Γ-centred 2x2x2 mesh (the 8
// TRIM points, listed), PADE, OUR basis valence_lowq_sr2_naq9.bsd (NAF-Q9-BASIS by scripts/bsd2cp2k; s and p only, so spherical =
// Cartesian), GTH-PADE q9 Na / q7 F, CUTOFF 1200 Ry.  Both start at a = 8.73 bohr.  CP2K's minimum: a 8.519577 bohr, E -72.0950276 Ha
// (6 steps); at 2400 Ry 8.519336 / -72.0950120, so CP2K's cutoff holds a to 2.4e-4 bohr.  Rocksalt is STABLE here: CP2K's free CELL_OPT
// from 0.002 fractional noise (`naf_q9_k222_free_cellopt.out`, 200 steps, its BFGS stalled) never went below -72.0948712.
namespace
{
//! The rocksalt NaF geometry run at cubic constant \a a0 (fcc primitive cell, Na q9 / F q7, PADE, k222, FREE density), with the
//! Fm-3m map over its group; runs \a body on the run's surface and the map.  The grids are STATED: those the default accuracy's raster
//! pilot settled on at a = 8.73 (2026-11-21: 3 audit rounds, the ball doubled 2 -> 4 -> 8, every term meeting 1e-5 Ha), as its "a deck
//! restates these grids as" line prints them.  The piloted relaxation reached the same a (8.519761) and E (-72.0950084) to 1e-9 Ha, but
//! its audit rounds (2.4 GB references) made it a 4 h run; stated, it is 2.4 h (OpenWork §4 "A GPW forces + stress point on NaF").
void WithNaF(double a0, const char* label,
             const std::function<void(qchem::PotentialEnergySurface&, const qchem::CellWyckoffMap&)>& body)
{
    using namespace qchem;
    const double h=0.5*a0;
    UnitCell cell(Matrix3D<double>(0.0,h,h,  h,0.0,h,  h,h,0.0));
    cell.AddAtom(11, rvec3_t(0.0,0.0,0.0));
    cell.AddAtom( 9, rvec3_t(0.5,0.5,0.5));
    SCFParams p; p.NMaxIter=100; p.MinΔρPerElectron=1e-7; p.MinΔE=std::nullopt; p.MinΔFD=std::nullopt; p.MinVirial=std::nullopt; p.MinFD=std::nullopt;
    RunOptions o; o.Nelec=16; o.species={{"Na",9},{"F",7}}; o.label=label; o.scf=p;
    o.model=Hamiltonian::Model::PADE;   // the GTH-PADE PPs' own functional, and CP2K's &PADE
    o.solid.cellStrainAllowance=0.05;   // the minimum is 2.4 % COMPRESSED from the cold start: the frozen sizes over-converge
    o.solid.xcMeshSymmetric=false;      // a FREE geometry run: one plain Becke mesh at every point
    o.accuracy.eps=0.0;                 // the pilot's grids, stated (above)
    o.solid.cutoffFactor=8.0;
    o.solid.tolerances.rasterX=o.solid.tolerances.xcRasterX=19.35139721329466;
    o.solid.tolerances.vlocRasterX=21.583068172270572;
    o.solid.tolerances.vlocEps=1e-6;
    o.xcMesh=qcMesh::BeckeXCParams(100, 2.0, 53);
    const auto setup=SolidSetup::Resolve(Lattice_3D(cell, ivec3_t(2,2,2)), {.file="valence_lowq_sr2_naq9.bsd"}, o);
    auto run=MakeSolidGeometryRun(setup, {SCFStage{p, SCFAccelerators::Type::Auto}}, nullptr, SolidSession::Carry::AO, /*relaxCell*/true);
    const PESGeometry ref=run.surface->Reference();
    const CellWyckoffMap map(ref, setup->lattice.GetSpaceGroup(), ElasticModel{}, 1e-4*a0);
    ASSERT_EQ(map.Dimension(), 1u) << "rocksalt's one free coordinate is the cubic constant";
    body(*run.surface, map);
}
} // namespace

//------------------------------------------------------------------------------------------------------------------------------------
// ROCKSALT NaF RELAXES ITS CELL ON GPW (doc/OpenWork.md §2 "Real relaxations end to end", the easy cases: controller plan 2026-11-20 (b)).
// Cell only: rocksalt has no free internal coordinate, so the map over the detected group (Fm-3m) has ONE coordinate, the cubic constant.
// FREE density (all 8 k-points, the plain Becke mesh), the group held by the map.  SEMICORE Na q9: with q1 Na the site has no restoring
// force, in CP2K too (OpenWork "A light-PP ion has NO restoring force"), so a q1 crystal is no stable structure to relax.  The oracle is
// CP2K CELL_OPT on the identical model (`~/Code/cp2k-runs/NaF/cellopt/naf_q9_k222_cellopt.{inp,out}`): the Γ-centred 2x2x2 mesh (the 8
// TRIM points, listed), PADE, OUR basis valence_lowq_sr2_naq9.bsd (NAF-Q9-BASIS by scripts/bsd2cp2k; s and p only, so spherical =
// Cartesian), GTH-PADE q9 Na / q7 F, CUTOFF 1200 Ry.  Both start at a = 8.73 bohr.  CP2K's minimum: a 8.519577 bohr, E -72.0950276 Ha
// (6 steps, 432 s on one rank); at 2400 Ry 8.519336 / -72.0950120, so CP2K's cutoff holds a to 2.4e-4 bohr.  Rocksalt is STABLE here:
// CP2K's free CELL_OPT from 0.002 fractional noise (`naf_q9_k222_free_cellopt.out`, 200 steps, its BFGS stalled) never went below
// -72.0948712.
// MEASURED 2026-11-21 (4 threads): 4 steps, 5 points, a 8.519761 (1.8e-4 bohr from CP2K), E -72.0950084 (19 uHa above CP2K's), 8659 s.
// So the cold run is in the PERIODIC tier.  Its twin, the ONE point below, is there too: one point costs 2328 s alone and timed out
// (3600 s) in the long tier beside two others (2026-11-22).  It returns to the long tier with OpenWork §4's Becke-force cost row.
TEST(GPW_NaF, k222_RelaxesToCP2KsCell_Long)
{
    using namespace qchem;
    WithNaF(8.73, "NaF relax", [](PotentialEnergySurface& pes, const CellWyckoffMap& map)
    {
        auto r=GeometryOptimizer(pes, map, OptimizerOptions::FromEps(1e-6)).Run();
        ASSERT_TRUE(bool(r)) << r.Error().details;
        report::FormatGuard keep(std::cout);
        std::cout<<"[relax] NaF: "<<r->Steps()<<" steps, "<<r->history.size()<<" evaluations"<<std::endl;
        for (const auto& it : r->history)
            std::cout<<"[relax]   "<<std::setw(3)<<it.index<<"  E "<<std::setprecision(10)<<std::fixed<<it.value<<std::defaultfloat
                     <<"  max|g| "<<std::setprecision(3)<<std::scientific<<it.maxGradient<<"  |s| "<<it.stepLength
                     <<std::defaultfloat<<(it.accepted ? "" : "  rejected")<<std::endl;
        const rmat3d_t& A=*r->geometry.cell;
        const double a=2.0*A(2,1);   // the fcc primitive vector (0, a/2, a/2) in column 1
        std::cout<<std::setprecision(10)<<"[relax] NaF a "<<a<<"  E "<<r->value<<"  (CP2K a 8.519577, E -72.0950276)"<<std::endl;
        EXPECT_NEAR(A(1,1), 0.0, 1e-8); EXPECT_NEAR(A(3,1), A(2,1), 1e-8);   // still fcc (the map holds it)
        EXPECT_NEAR(a, 8.519577, 5e-3);
        EXPECT_NEAR(r->value, -72.0950276, 1e-4);
    });
}

// THE TWIN: CP2K's minimum is a stationary point of OUR surface -- one point, the map gradient at CP2K's cubic constant.  The cold start
// (a 8.73) has max|g| 7.5e-3; the relaxation above ends 1.8e-4 bohr from this point.
TEST(GPW_NaF, k222_CP2KsCellIsStationary_Long)
{
    using namespace qchem;
    WithNaF(8.519577, "NaF twin", [](PotentialEnergySurface& pes, const CellWyckoffMap& map)
    {
        const rvec_t q=map.Start();
        const PESGeometry g=map.At(q);
        auto point=pes.ConvergeAt(g.cell, g.R, map.PointRequest());
        ASSERT_TRUE(bool(point)) << point.Error().details;
        const rvec_t grad=map.Project(q, map.Gradient(q, *point));
        report::FormatGuard keep(std::cout);
        std::cout<<std::setprecision(10)<<"[stationary] NaF at CP2K's cell: E "<<point->energy<<" (CP2K -72.0950276)  g "
                 <<std::setprecision(3)<<std::scientific<<grad[0]<<std::defaultfloat<<std::endl;
        EXPECT_LT(std::fabs(grad[0]), 1e-4);
        EXPECT_NEAR(point->energy, -72.0950276, 1e-4);
    });
}

// THE PERIODIC GGA ON NaF (GGA stage 3): the free Becke delta route, Pulay on D (the default mixer: a GGA refuses Kerker), SR2, Gamma, PBE
// (gth-pbe-valence), the XC mesh UNSTATED: the GGA route chooses Becke and the pilot sizes it from eps (controller ruling 2026-11-21, pin 31).
// The oracle is CP2K naf_gpw_sr2_diag.inp at PBE, converged in CUTOFF (~/Code/cp2k-runs/NaF/naf_sr2_pbe*): 320 / 600 / 1000 / 1600 Ry =
// -24.453123 / -24.452611 / -24.452574 / -24.4525727; libxc's PBE at 1000 Ry -24.4525716, XC cutoffs 1e-14 -24.4525740.
// WHY THE MESH IS NOT STATED (measured 2026-11-21, ours at Becke nR 150, LDA / PBE): L47 -24.4310540 / -24.4525424, L65 -24.4310536 /
// -24.4525718 -- the GGA moved 29 uHa between two angular rungs that moved LDA 0.4 uHa, so the energy-calibrated default (40, 29) cannot
// serve a GGA; radial 150 -> 200 at L47 moved PBE 1.2 uHa.
TEST(GPW_NaF, Γ_Becke_PBE_CP2K_Long)
{
    const Material naf=qchem::Materials::Get("NaF_rocksalt");
    const Lattice_3D lat=LatticeOf(naf);
    RunOptions o=NaFOptions(naf, "NaF Gamma Becke PBE");
    o.model=qchem::Hamiltonian::Model::PBE; o.pp=std::string("gth-pbe-valence");
    SCFParams par=NaFGates(); par.KerkerG0=0.0; par.StartingRelaxRo=0.45;
    o.scf=par;
    auto calc=qchem::SolidCalculation::Prepare(SolidSetup::Resolve(lat, MakeBasisNaFSR2(*naf.cell), o)).TakeValue(); (void)calc->Converge();
    auto R=calc->Result();
    ASSERT_TRUE(R) << Why(R);
    EXPECT_NEAR(R->TotalCharge(), 8.0, 1e-6);
    EXPECT_NEAR(R->Energy(), -24.4525727, 2e-5) << "PBE on the mesh the run sized from eps (1e-5 x 2 atoms) vs CP2K at 1600 Ry";
}
