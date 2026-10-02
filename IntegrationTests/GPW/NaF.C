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

#include "gtest/gtest.h"
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
#include <iomanip>      // setprecision (the order-parameter trajectory line)

import qchem.Structure;                          // Molecule, Atom
import qchem.UnitCell;                           // UnitCell, FCCUnitCell
import qchem.Matrix3D;                           // Matrix3D<double> (the rhombohedral AFM-II cell matrix)
import qchem.Lattice_3D;                         // Lattice_3D
import qchem.BasisSet;                           // Complex_BS, Real_BS
import qchem.BasisSet.Orbital_1E_IBS;            // Complex_OIBS (the overlap-spectrum diagnostic)
import qchem.Blaze;                              // blazem::eigen, blaze::min/max (overlap spectrum)
import qchem.BasisSet.PlaneWave.PlaneWave_IBS;   // PlaneWave_IBS (the seed's CD fit basis)
import qchem.BasisSet.Lattice.BasisSet;       // GPWFactory (the GPW basis container)
import qchem.BasisSet.Gaussian.Point.Factory;          // Gaussian::Factory, BasisSetData/Engine/Angular
import qchem.BasisSet.Gaussian.Lattice.SphericalLatticeView;  // MakeSphericalLatticeView (GPW_SPHERICAL=1)
import qchem.Hamiltonian.Factory;                 // the PUBLIC solid front door (Step 4): cHamiltonian* Factory(...)
import qchem.Outcome;                           // Outcome<Converged,SCFFailure> -- the facade's result
import qchem.RunPolicy;                         // ReresolveRunPolicy() -- the declared-deviation A/B hatch (N5)
import qchem.SolidCalculation;                    // the NAMED periodic facade (Step 4 3/3)
import qchem.Tests.GPW_Harness;                   // THE HARNESS (IntegrationTests/GPW/Harness.C): Materials cells, gates, recipes, the XC probes
import qchem.Materials;                           // Materials::Get -- the cells come from src/Structure/Data/materials.json (row MD)
import qchem.Hamiltonian.Internal.Hamiltonians;  // Ham_PW_DFT direct ctors (the bespoke probes below still use them)
import qchem.Hamiltonian.Internal.PWTerms;        // ReportGridCharge(); Vxc_Quadrature + the two DensitySampler strategies
import qchem.ChargeDensity.DensitySampler;  // the XC sampling engine (its own module
                                                  // since 2026-09-08; .Internal. modules are
                                                  // never re-exported, so name it directly)
import qchem.BasisSet.DeltaFit_IBS;              // DeltaFit_IBS -- the delta basis the singles strategy runs on
import qchem.BasisSet.G_FieldEvaluator;           // G_RasterTransform -- the uniform probe's own point count
import qchem.Mesh.Angular;                        // MakeAngular (the rotated-Lebedev bond-angle probe)
import qchem.Hamiltonian.Internal.ExFunctional;   // ExFunctional (the LDA functional face the XC terms hold)
import qchem.Hamiltonian.Internal.SlaterExchange; // SlaterExchange (Dirac exchange, for the Becke XC gate)
import qchem.Hamiltonian.Internal.VWN_Correlation;// VWN_Correlation (VWN5, for the Becke XC gate)
import qchem.Mesh;                                // qcMesh::MeshParams / UnitCellKind (the Becke XC quadrature)
import qchem.Mesh.XCPolicy;                       // BeckeXCParams / ResolveXCMesh / XCMeshSharpness (the grid policy)
import qchem.BasisSet.Gaussian.Lattice.LatticeSum1E;      // Gaussian::LatticeSum1E::MaxExponent (alpha_max, for the selector)
import qchem.Pseudopotential.GTH_Potentials;      // GetGTH -> HGH local PP (alpha_pp = 1/2r_loc^2, for the selector)
import qchem.PeriodicTable;                       // thePeriodicTable().GetZ (element symbol -> Z)
import qchem.SCFIterator;                        // cSCFIterator, SCFParams
import qchem.SCFParams;                          // SCFParams
import qchem.ElectronConfiguration.Crystal;      // Crystal_EC (single-k Bloch occupation)
import qchem.ChargeDensity.Seed;                 // SeedStrategy
import qchem.SCFAccelerator.Factory;              // the PUBLIC complex accelerator door (Step 4)
import qchem.SCFAccelerator.Internal.SCFAcceleratorDIIS; // SCFAcceleratorDIIS (scalar-agnostic manager)
import qchem.SCFAccelerator.Internal.SCFAcceleratorGDM;  // SCFAcceleratorGDM (scalar-agnostic manager)
import qchem.SCFAccelerator.Internal.SCFAcceleratorLadder; // SCFAcceleratorLadder (DIIS -> GDM chain)
import qchem.SCFAccelerator.Internal.SCFIrrepAcceleratorNull; // SCFAcceleratorNull (NaF: pure damped Kerker)
import qchem.WaveFunction;                       // cWaveFunction (the converged state)
import qchem.Energy;                             // EnergyBreakdown
import qchem.Symmetry.Irrep;                     // Irrep
import qchem.Reporting;                          // report:: -- bracket the GPW run so grids/basis sections land
import qchem.Symmetry.Spin;                      // Spin
import qchem.Symmetry.Factory;                   // BlochFactory (build a k-block with a fractional MP shift)
import qchem.LASolver;                           // qchem::Ortho (Cholesky | Eigen | SVD -- basis orthogonalisation)
import qchem.BasisSet.Gaussian.Lattice.GPW_IBS;         // GPW_IBS (build a concrete block for the collocation diagnostic)
import qchem.BasisSet.Gaussian.Lattice.GPW_Evaluator;  // GPW_Evaluator (Overlap3CTensor -- the collocation tensor)
import qchem.BasisSet.GMap;              // Projector3<dcmplx> (the collocation weight tensor); SymmetryDefects (§3 diagnostic)
import qchem.ChargeDensity.FourierDensity;        // FourierDensity (ρ̃ for the §3 order-parameter diagnostic)
import qchem.CompositeCD;                         // tComposite_CD (the polarized density = one composite over Up+Down blocks, V1.37)
import qchem.ChargeDensity.Factory;
import qchem.ChargeDensity.SeedCD;              // PolarizedSeedCD (the raw spin-SAD seed, for the sublattice gate)               // IrrepCD_Factory/PolarizedCD_Factory (fixed-density probe)
import qchem.Pseudopotential.GTH_Potentials;      // GetGTH, GTH_PP (the PP model, for the matrix-trace probe)
import qchem.Calculation;                        // qchem::Calculation, CalcOptions (finite reference)
import qchem.AtomCalculation;                    // AtomCalculation, AtomType, BasisSetAccuracy (Slater/High pseudo-atom ref)
import qchem.Types;

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
// Explicit densityEcut= 40 = SUB-FLOOR: warns, and BallOnly aliases there (-43 mHa)
TEST(GPW_NaF, Γ_Imp_Anchor)   // RE-ENABLED 2026-09-15: 15 s, converged in 23 iterations -- the GPW x NaF row's first standing anchor
{
    // THE NaF ORACLE ANCHOR (doc/GPWPlan.md: 0.10-0.19 mHa vs CP2K on this basis, tight-eps + converged density).
    // The committed production recipe (NaFOptions/NaFGates), SR2 basis, Γ -- the arm doc/Benchmark.md times and
    // the CP2K decks (naf_gpw_sr2_diag.inp / naf_gpw_sr_tight.inp, no &KPOINTS) compute.  The k-mesh, span,
    // ecut, ladder, alpha, smearing and penalty knobs this test used to read from NAF_* env vars are the
    // campaign INSTRUMENT's and go with it to the probe binary (doc/TestSuitePlan.md §8).
    const Material naf=qchem::Materials::Get("NaF_rocksalt");
    const Lattice_3D lat=LatticeOf(naf);
    SolidCalcOptions o=NaFOptions(naf, "NaF GPW Gamma");
    o.imposeSymmetry=true;   // V1.30: was the DEFAULT; now stated, because an imposition you did not ask for is invisible in the result
    o.xcMesh          = qcMesh::BeckeXCParams(20,2,24);
    o.xcMesh.cellKind = qcMesh::UnitCellKind::Becke;
    SCFParams par=NaFGates(); par.StartingRelaxRo=0.25; par.MergeTol=1e-4;
    par.Verbose=(bool)std::getenv("GPW_VERBOSE");
    EnvOverrides(o, par);                            // the Benchmark row (scripts/retake5a) drives it
    GpwReport report("NaF "+o.label, par.Verbose);
    qchem::SolidCalculation calc(lat, MakeBasisNaFSR2(*naf.cell), o, par);
    auto R=calc.Result();
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
    bool        ok=false;
    std::string why;
};
SolidCalcOptions NaFAnchorOptions(const Material& naf, const std::string& label)
{
    SolidCalcOptions o=NaFOptions(naf, label);
    o.imposeSymmetry=true;
    o.xcMesh          = qcMesh::BeckeXCParams(20,2,24);
    o.xcMesh.cellKind = qcMesh::UnitCellKind::Becke;
    return o;
}
SCFParams NaFAnchorParams() { SCFParams par=NaFGates(); par.StartingRelaxRo=0.25; par.MergeTol=1e-4; par.Verbose=(bool)std::getenv("GPW_VERBOSE"); return par; }

const NaFRestartFixture& NaFRestart()
{
    static const NaFRestartFixture f=[]
    {
        NaFRestartFixture r;
        const Material naf=qchem::Materials::Get("NaF_rocksalt");
        const Lattice_3D lat=LatticeOf(naf);
        const auto dir=std::filesystem::path(testing::TempDir());
        r.coarsePath=(dir/"naf_restart_coarse.h5").string();
        r.finePath  =(dir/"naf_restart_fine.h5").string();
        {   // the CHEAP coarse stage: Ecut=40 is SUB-FLOOR (below C*alpha_max=80), so it needs the alias-free raster
            SolidCalcOptions o=NaFAnchorOptions(naf, "NaF restart coarse");
            o.densityEcut=40.0;
            o.raster=BasisSet::PlaneWave::RasterPolicy::AliasFree;
            o.saveStateTo=r.coarsePath;
            qchem::SolidCalculation calc(lat, MakeBasisNaFSR2(*naf.cell), o, NaFAnchorParams());
            auto R=calc.Result();
            if (!R) { r.why="coarse: "+Why(R); return r; }
            r.coarseEnergy=R->Energy();
        }
        {   // the COLD fine reference (what Γ_Imp_Anchor runs), also saved for the exact-resume claim
            SolidCalcOptions o=NaFAnchorOptions(naf, "NaF restart cold fine");
            o.saveStateTo=r.finePath;
            qchem::SolidCalculation calc(lat, MakeBasisNaFSR2(*naf.cell), o, NaFAnchorParams());
            auto R=calc.Result();
            if (!R) { r.why="cold fine: "+Why(R); return r; }
            r.coldEnergy=R->Energy(); r.coldIterations=R->IterationCount();
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
TEST(GPW_NaF, Γ_Imp_eqColdStart)
{
    const NaFRestartFixture& f=NaFRestart();
    ASSERT_TRUE(f.ok) << f.why;
    EXPECT_NEAR(f.coarseEnergy, -24.4357, 0.01) << "the coarse seed stage must itself be the physical fixed point";
    const Material naf=qchem::Materials::Get("NaF_rocksalt");
    const Lattice_3D lat=LatticeOf(naf);
    SolidCalcOptions o=NaFAnchorOptions(naf, "NaF restart fine");          // auto densityEcut = the production grid
    std::vector<double> E;
    o.onIteration=[&E](const qchem::SCFIterator::SCFProgress& p){ E.push_back(p.energy); };
    auto c=qchem::SolidCalculation::Restart(f.coarsePath, lat, MakeBasisNaFSR2(*naf.cell), o, NaFAnchorParams());
    ASSERT_TRUE(c) << (c ? std::string() : c.Error().details);
    auto R=(*c)->Result();
    ASSERT_TRUE(R) << Why(R);
    ASSERT_FALSE(E.empty());
    EXPECT_NEAR(R->TotalCharge(), 8.0, 1e-6);
    EXPECT_NEAR(R->Energy(), f.coldEnergy, 1e-4) << "the restarted run reached a different state than the cold fine run";
    EXPECT_NEAR(R->Energy(), -24.4304, 0.01) << "and it is the banked NaF anchor";
    EXPECT_LT(std::abs(E.front()-f.coldEnergy), 0.05) << "the first fine iterate should already be near the answer (coarse seed)";
    EXPECT_LT(R->IterationCount(), f.coldIterations) << "the coarse seed must SAVE iterations on the fine grid (measured 12 vs 21)";
}

// CLAIM: the SAME grid resumes EXACTLY -- Restart from the converged fine state reproduces its energy and needs
// (almost) no iterations.  The control for the claim above: if this fails the file round-trip is lossy on a hard
// material, and the grid-continuation result means nothing.
TEST(GPW_NaF, Γ_Imp_eqExactResume)
{
    const NaFRestartFixture& f=NaFRestart();
    ASSERT_TRUE(f.ok) << f.why;
    const Material naf=qchem::Materials::Get("NaF_rocksalt");
    const Lattice_3D lat=LatticeOf(naf);
    SolidCalcOptions o=NaFAnchorOptions(naf, "NaF restart exact");
    o.accelerator=qchem::SCFAccelerators::Type::GDM;   // a converged state wants GDM, not the DIIS+Kerker rung (see below)
    std::vector<double> E;
    o.onIteration=[&E](const qchem::SCFIterator::SCFProgress& p){ E.push_back(p.energy); };
    auto c=qchem::SolidCalculation::Restart(f.finePath, lat, MakeBasisNaFSR2(*naf.cell), o, NaFAnchorParams());
    ASSERT_TRUE(c) << (c ? std::string() : c.Error().details);
    auto R=(*c)->Result();
    ASSERT_TRUE(R) << Why(R);
    ASSERT_FALSE(E.empty());
    // GDM, NOT the Ladder (user: restarts want GDM, or the Ladder's GDM rung).  MEASURED 2026-10-02: from the same file
    // the Ladder's DIIS rung + Kerker(alpha=.25) STEPS the density before the first energy is read -- first iterate
    // 1.3e-4 Ha off the saved state, 4 DIIS iterations until |dE/E|<1e-6 hands off to GDM, 9 in all -- whereas GDM
    // starts at the answer (4e-10) and converges in 2.  The saved state was never the problem.
    EXPECT_NEAR(E.front(), f.coldEnergy, 1e-7) << "the first iterate of an exact resume is the saved state";
    EXPECT_LE(R->IterationCount(), 3u) << "an exact resume under GDM starts converged";
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
    // BallOnly): the SCF only supplies the density; the comparison itself carries the reference-grade work.
    // MEASURED ATTRIBUTION (2026-07-30): vs the SAME Becke B40 matrix, the uniform reference gave
    //   BallOnly  Ecut=160: max|U-B|=1.55e-1      (the production raster's sharp-F element error)
    //   AliasFree Ecut=160: max|U-B|=1.89e-2      (8x better -- most of the gap was BallOnly's)
    // with dExc tiny throughout (1.6e-4 / 1.3e-5) -- the discrepant elements carry no rho weight.  So on a
    // sharp-F system the UNIFORM raster is not element-converged at practical Ecut, which is this grid's
    // reason to exist; the gate here is Becke INTERNAL convergence (B40 vs a 2x-refined B80) plus the
    // energy-level agreement with the uniform route.
    SolidCalcOptions o=NaFOptions(naf, "NaF Becke gate");
    o.imposeSymmetry=true;   // V1.30: was the DEFAULT; now stated, because an imposition you did not ask for is invisible in the result
    GpwReport report("NaF "+o.label, false);
    qchem::SolidCalculation calc(lat, MakeBasisNaFSR2(*naf.cell), o, NaFGates());
    ASSERT_TRUE(calc.Result()) << Why(calc.Result());
    EXPECT_NEAR(calc.Result()->TotalCharge(), 8.0, 1e-6);
    const GpwHandles h=Handles(calc);

    auto st=lat.GetStructure();
    XCProbe U  =UniformXCProbe(h, st);
    XCProbe B40=BeckeXCProbe(h, st, "B40", qcMesh::BeckeXCParams());
    XCProbe B80=BeckeXCProbe(h, st, "B80", qcMesh::BeckeXCParams(/*nRadial*/80, /*mhlAlpha*/1.0, /*L*/29));
    EXPECT_NEAR(B40.Exc, U.Exc, 2e-3);               // energy-level agreement with the uniform route
    EXPECT_NEAR(B40.rhoLost, 0.0, 5e-3);
    DiffXC(U,B40);                                    // report-only: the uniform raster's element error
    EXPECT_LT(DiffXC(B40,B80), 2e-3) << "Becke V_xc not internally converged on the sharp-F system";
}


// (The pre-facade grid-continuation test, DISABLED_Γ_GridContinuation -- raw SolidSCFIterator ctors, ~15 GC_* env
// knobs -- was DELETED 2026-10-02 (D-NAFGRID): its basin-avoidance premise (the direct fine run diving to -40 Ha) no
// longer reproduces, and its incremental-convergence claim lives in Γ_Imp_eqColdStart / Γ_Imp_eqExactResume above.
// History: doc/Records/OpenWork_History*.md, doc/OldPlans/GPWPlan.md §0e.)
