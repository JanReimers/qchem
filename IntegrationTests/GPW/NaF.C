// File: IntegrationTests/GPW/NaF.C  NaF rocksalt -- the sharp-F ionic multi-species cell (a=8.73, SR2 basis); the best CP2K-validated material (0.2 mHa).
//
// THE GRID (doc/TestSuitePlan.md §3): TEST(<Basis>_<Material>, <k>_[Grid]_[Fit]_[Sym]_[Spin]_[Occ]_[Machinery]_[Ansatz]_[Seed]_<Claim>)
// -- axis tokens in fixed order, the facade's defaults elided (k always named); the claim is CP2K (oracle anchor),
// Anchor (did-E-move), eq<Token> (a ONE-axis twin: this point equals the same point with that axis moved) or a
// property verb.  Tests are laid out in axis order.  `scripts/testgrid` renders the coverage table from these names.
//
//   GPW_NaF.Γ_Imp_Anchor
//   GPW_NaF.Γ_Imp_Becke_eqUni
//   GPW_NaF.DISABLED_Γ_GridContinuation

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
import qchem.Materials;                           // Materials::Get -- the cells come from src/Calculation/Data/materials.json (row MD)
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
    GpwReport report("NaF "+o.label, par.Verbose);
    qchem::SolidCalculation calc(lat, MakeBasisNaFSR2(*naf.cell), o, par);
    auto R=calc.Result();
    ASSERT_TRUE(R) << Why(R);
    EXPECT_NEAR(R->TotalCharge(), 8.0, 1e-6);           // 1 (Na) + 7 (F) valence electrons, conserved
    EXPECT_NEAR(R->Energy(), -24.4304, 0.01);           // did-E-move anchor at Γ (CP2K agrees to 0.2 mHa, doc/GPWPlan.md)
}



// The SHARP-FIELD leg of the gate (the plan names DISABLED_NaFRocksaltGamma as the stress case: the F-
// anion makes sharp peaks in rho and V_xc, and its diffuse basis is what the Becke grid exists for).
// DISABLED like the parent NaF anchor -- it is a long run (the full NaF convergence recipe at a
// matrix-grade densityEcut=160 reference, plus two Becke term evaluations); run it by hand with
// --gtest_also_run_disabled_tests when touching the XC quadrature.
TEST(GPW_NaF, Γ_Imp_Becke_eqUni)   // RE-ENABLED 2026-09-15: 24 s, Becke internally converged on the sharp-F system
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


// (4b) NaF GRID-CONTINUATION SEEDING (doc/GPWPlan §0e, step 1) -- the PRODUCTION-GRID fix.
// (Pins + story re-derived 2026-07-23 on the post-analytic-short/kappa/5-smooth landscape; the original
// -27.76/-27.93 anchors were the RETRACTED aliasing-era values -- doc/GPWPlan.md TRAPS #2.  The clean NaF
// SR2 truth is CP2K -24.4312 at tight eps, Ecut=160-class grids.)
//
// THE PROBLEM.  The direct production-grid NaF run FALLS INTO the unphysical XC-collapse basin: from the
// ionic seed, the Kerker priming descent goes straight into E~-40 (mid-slosh D loads the sharpest F pairs
// beyond the grid calibration -> the collocated rho aliases spiky/locally-negative -> E_xc is legitimately
// huge-negative WITHIN the discretization, a self-consistent garbage fixed point), and Pulay engaging on
// that garbage state thrashes to +54.  MOM+Pulay are NECESSARY but NOT SUFFICIENT: the basin is a property
// of the map, reachable by the descent -- not an occupation swap (MOM) nor a mixing wobble (Pulay).
// HISTORY: on the BALL-XC map the basin was real at every sub-C=8 grid (the 0.5(f1) sweep hit it from a
// SEEDED start at Ecut=160: negCharge -91, Exc -109).  The 0.5(f2) raw-XC feed REMOVED it (rho_DM >= 0 by
// construction -- negCharge == 0 at C=8/4/3 in the acceptance sweep); this test retains the ionic-seed A/B
// as the historical repro knob.
//
// THE FIX (this test).  Never ENTER the basin: converge the CHEAP coarse grid (Ecut=40), then SEED the fine
// grid with that converged density so the fine SCF STARTS in the physical basin.  The orbital (SR2 Gaussian)
// basis is IDENTICAL at both cutoffs -- only the density COLLOCATION grid differs -- so the converged coarse
// density matrix transfers directly (no re-projection) via the explicit-seed cSCFIterator ctor, which
// collocates the seed's iteration-0 Hartree/XC on the REQUESTED fine fit grid (the fit-grid seam is honest
// since 2026-07-20 -- GPW_IBS builds the tensor over the requested fit basis's grid).  Init immediately
// re-diagonalizes on the fine basis and every subsequent iteration runs the fine grid, starting in (and
// staying in) the physical basin.
//
// MOM ACROSS THE GRID CHANGE (doc/GPWPlan 0h): transferring the coarse WF's occupied subspace as a fixed
// MOM reference (AdoptMOMReference) pinned AN EXCITED STATE across the discretization change (measured
// 2026-07-23: -23.680, +0.75 Ha).  The 0h GUARD (persistent-hole detection -> release + re-capture) now
// makes BOTH recipes land the ground state: pure aufbau (the default: GC_SEED_MOM=0, GC_FINE_MOM_START=
// 9999; 22 iters) and the transfer path (GC_SEED_MOM=1 GC_FINE_MOM_START=1: VERIFIED 2026-07-23, guard
// fires once on the coarse stage's own capture-at-10 reference, fine converges 16 iters to -24.43252 --
// identical to the aufbau pin to 8 decimals).  The guard also exposed that the COARSE stage's endpoint had
// itself been MOM-pinned +0.75 high in every earlier measurement (see the coarse-pin note below).
//
// GATE: the fine grid must reach the raw-XC aufbau ground state -24.4325 (1.3 mHa from CP2K's -24.4312,
// itself an Ecut=160-class number), NOT the -40 basin.  DISABLED (two full NaF SCFs, ~5 min).  Env knobs
// (GC_*) tune each stage without recompiling.  Verify basin-avoidance is REAL by A/B: with GC_SEED=0 the
// fine stage falls back to the ionic seed and must dive into the basin (the direct-run failure this test
// fixes; the energy gates then fail by design).
TEST(GPW_NaF, DISABLED_Γ_GridContinuation)
{
    using namespace qchem::Hamiltonian;
    namespace L3=BasisSet::Lattice;
    auto envd=[](const char* n, double d){ const char* s=std::getenv(n); return s ? std::atof(s) : d; };

    const double a=8.73;
    FCCUnitCell cell(a);
    cell.AddAtom(11, {0,0,0});          // Na (Zion=1)
    cell.AddAtom(9,  {0.5,0.5,0.5});    // F  (Zion=7)
    Lattice_3D lat(cell, ivec3_t(1,1,1));
    auto st = lat.GetStructure();       // held for both stages (the ctors' non-owning structure view)
    // GC_BASIS selects the orbital basis (default SR2 = the well-conditioned regression config; "SR" = the
    // FULL short-range basis, lambda_min~1e-6 at complete enumeration -- the sec-1 rank-reduction campaign's
    // probe target, oracle CP2K -27.93128 on VALENCE-LOWQ-SR).
    const char* gcb=std::getenv("GC_BASIS");
    const BasisSetData basis = (gcb && std::string(gcb)=="SR") ? BasisSetData::VALENCE_LOWQ_SR
                                                               : BasisSetData::VALENCE_LOWQ_SR2;
    auto mol = std::shared_ptr<const Real_BS>(BasisSet::Gaussian::Factory(
        basis, &cell, BasisSet::Gaussian::Engine::MnD, BasisSet::Gaussian::Angular::Cartesian));

    // The converged Ecut=40 recipe (DISABLED_NaFRocksaltGamma): pure damped Kerker (NO DIIS), exit on E-flat,
    // delayed-IMOM MOM + Kerker-preconditioned Pulay.  Tunable per stage (near the fixed point the fine stage
    // can capture MOM + engage Pulay earlier, since it does not need the ~10/35-iter descent the coarse one does).
    auto makePar=[&](size_t nmax, int momStart, int pulayStart)
    {
        SCFParams par; par.NMaxIter=nmax; par.MinΔρ=1e30; par.MinΔE=1e-8; par.MinΔFD=1e30; par.MinVirial=1e30;
        par.MinFD=1e30; par.StartingRelaxRo=envd("GC_ALPHA",0.025); par.MergeTol=1e-4; par.Verbose=true;
        par.KerkerG0=envd("GC_KERKER_G0",1.0);
        par.UseMOM=true; par.MOMStartIter=momStart; par.PulayDepth=(int)envd("GC_PULAY",6); par.PulayStart=pulayStart;
        return par;
    };

    // RSS breadcrumb (the full-SR allocation-bomb bisect, 2026-07-22): prints resident MB per ctor phase.
    auto rss=[](const char* tag)
    {
        std::ifstream f("/proc/self/statm"); size_t vmpg=0, rspg=0; f>>vmpg>>rspg;
        std::cerr<<"[rss] "<<tag<<": "<<(rspg*4096/1048576)<<" MB"<<std::endl;
    };
    // ---- STAGE 1: converge on the CHEAP coarse density grid (Ecut=40 -> the physical fixed point). ----
    // Every coarse-stage object is a unique_ptr so the WHOLE stage can be torn down mid-test (below) the
    // moment the fine stage has consumed it -- doc/GPWPlan.md 0.5(b).
    rss("pre-basis");
    // The coarse SEED stage runs Ecut=40 -- SUB-FLOOR (below C*alpha_max=80), where BallOnly aliases
    // (-43 mHa); pin it to the exact-quadrature raster so the seed is the honest -24.4357 fixed point.
    std::unique_ptr<Complex_BS> bsC(L3::GPWFactory(lat, mol,
        L3::GPWParams{.densityEcut=envd("GC_COARSE_ECUT",40.0), .raster=BasisSet::Gaussian::RasterPolicy::AliasFree}));
    rss("basis");
    auto ecC=std::make_unique<Crystal_EC>(bsC->GetIrreps(Spin::None), 8);
    rss("EC");
    cHamiltonian* hamC=new Ham_PW_DFT(st, bsC.get(), {{"Na",1},{"F",7}}, "LDA");
    std::unique_ptr<cHamiltonian> hamCOwner(hamC);   // R2.22: the iterator borrows; this scope owns
    rss("Ham");
    auto* accC=new qchem::SCFAccelerators::SCFAcceleratorNull();   // no DIIS (the CP2K recipe)
    std::unique_ptr<qchem::SCFAccelerators::SCFAccelerator> accCOwner(accC);   // R2.22: the iterator borrows; this scope owns
    auto scfC=std::make_unique<qchem::SCFIterator::SolidSCFIterator>(bsC.get(), ecC.get(), hamC, accC,
                                          qchem::ChargeDensity::SeedStrategy::IonicSAD, st.get(),
                                          qchem::Cholesky, 0.0);
    rss("SCFctor");
    qchem::ChargeDensity::ReportGridCharge()=(bool)std::getenv("GPW_GRIDCHARGE");   // step-2 probe: coarse-grid rho stats to compare vs fine
    scfC->Iterate(makePar((size_t)envd("GC_COARSE_NMAX",200), 10, 35));
    qchem::ChargeDensity::ReportGridCharge()=false;
    auto Ecoarse=scfC->GetEnergy();
    std::cout << "[NaF grid-cont COARSE] Ecut=40 iters="<<scfC->GetIterationCount()
              << " Etot="<<Ecoarse.GetTotalEnergy() << std::endl;
    // The Ecut=40 fixed point on the RAW-XC landscape WITH the 0h MOM guard: -24.4357, E-flat converged in
    // ~43 iters -- only 3.2 mHa from the fine (Ecut=320) -24.4325: under raw XC the Ecut=40 grid is nearly
    // converged.  PIN HISTORY (each anchor exposed by the next fix): -27.76 = the RETRACTED aliasing era;
    // -23.69 (ball, 515 iters) and -23.68 (raw, 45 iters) = a MOM-PINNED EXCITED STATE the 0h guard caught
    // (the capture-at-fill-10 reference grabbed a non-aufbau configuration; persistent ~3 mHa hole ->
    // release -> aufbau recovery).  The "coarse underbinds by 0.74 Ha" story was that excited state's
    // artifact, not grid error.
    EXPECT_NEAR(Ecoarse.GetTotalEnergy(), -24.4357, 0.01);   // seed-quality anchor (did-E-move)

    // Grab the converged coarse density (OWNED; consumed by the fine ctor's Init).  bsC stays alive until
    // after the fine ctor, so the density's coarse-block pointer stays valid for the one iteration-0 read.
    auto* seedCD = scfC->GetWaveFunction()->GetChargeDensity().release();   // consumed by the fine ctor

    // ---- STAGE 2: seed the PRODUCTION fine grid (auto Ecut=8*alpha_max=320) with the converged coarse density. ----
    std::unique_ptr<Complex_BS> bsF(L3::GPWFactory(lat, mol, /*densityEcut*/envd("GC_FINE_ECUT",-1.0)));  // <0 AUTO=320
    Crystal_EC ecF(bsF->GetIrreps(Spin::None), 8);
    // R2.22: the iterator borrows these; this scope owns them, and they outlive scfF below.
    std::unique_ptr<cHamiltonian> hamF(new Ham_PW_DFT(st, bsF.get(), {{"Na",1},{"F",7}}, "LDA"));
    std::unique_ptr<qchem::SCFAccelerators::SCFAccelerator> accF(new qchem::SCFAccelerators::SCFAcceleratorNull());
    qchem::ChargeDensity::ReportGridCharge()=(bool)std::getenv("GPW_GRIDCHARGE");
    qchem::SCFIterator::ReportBandGap()=true;
    // GC_SEED=0 A/Bs the fix OFF (ionic seed) -> the fine stage must dive into the -39 basin (the failure this
    // test fixes); default ON = the converged-coarse-density explicit seed.
    const bool useSeed = envd("GC_SEED",1.0)!=0.0;
    std::unique_ptr<qchem::SCFIterator::cSCFIterator> scfF;
    if (useSeed)
    {
        scfF.reset(new qchem::SCFIterator::SolidSCFIterator(bsF.get(), &ecF, hamF.get(), accF.get(), seedCD, st.get(),
                                                        qchem::Cholesky, 0.0));   // explicit-seed ctor (consumes seedCD)
        // MOM transfer across the grid change is OFF by default: AdoptMOMReference across a discretization
        // change PINS AN EXCITED STATE (doc/GPWPlan 0h; measured 2026-07-23: -23.680 vs the -24.434 aufbau
        // ground state, +0.754 Ha).  With the density seed holding the run in the physical basin, the pure
        // aufbau fill converges cleanly to the ground state.  GC_SEED_MOM=1 (with GC_FINE_MOM_START=1)
        // re-enables the transfer path for the 0h MOM-guard work.
        if (envd("GC_SEED_MOM",0.0)!=0.0)
            scfF->AdoptMOMReference(*scfC->GetWaveFunction());
    }
    else
    {
        delete seedCD;   // A/B control: discard the coarse density, fall back to the ionic seed (dives to -39)
        scfF.reset(new qchem::SCFIterator::SolidSCFIterator(bsF.get(), &ecF, hamF.get(), accF.get(),
                                                        qchem::ChargeDensity::SeedStrategy::IonicSAD, st.get(),
                                                        qchem::Cholesky, 0.0));
    }
    // The coarse stage is DONE (seed consumed by the fine ctor's Init, MOM reference copied out) -- tear it
    // down IN ORDER (iterator -> EC -> basis) BEFORE the fine iterations (doc/GPWPlan.md 0.5(b)).  The
    // coarse ~GPW_Evaluator hands its ladder's stream caches back to the global budget, and the fine shape
    // (built STARVED during the handoff, while the coarse caches were still resident) rebuilds into the
    // refunded budget at its next EnsureStreams (the self-heal).  Without this the fine stage runs at ~0%
    // stream coverage, re-evaluating billions of points per iteration (the 8.45-h full-SR run).
    scfC.reset(); ecC.reset(); bsC.reset();
    rss("coarse stage freed");
    scfF->Iterate(makePar((size_t)envd("GC_FINE_NMAX",100),
                          (int)envd("GC_FINE_MOM_START",9999), (int)envd("GC_FINE_PULAY_START",12)));
    qchem::ChargeDensity::ReportGridCharge()=false;
    qchem::SCFIterator::ReportBandGap()=false;

    auto Efine=scfF->GetEnergy();
    auto cd=scfF->GetWaveFunction()->GetChargeDensity(); double charge=cd->GetTotalCharge();
    std::cout << "[NaF grid-cont FINE] auto-Ecut iters="<<scfF->GetIterationCount()<<" charge="<<charge
              << " Etot="<<Efine.GetTotalEnergy()
              << " (Ekin="<<Efine["Kinetic"]<<" Een="<<Efine["Een"]<<" Eee="<<Efine["Eee"]<<" Exc="<<Efine["Exc"]
              << " Enn="<<Efine["Enn"]<<" E_alphaZ="<<Efine["E_alphaZ"]<<")" << std::endl;
    EXPECT_NEAR(charge, 8.0, 1e-6);     // 1 (Na) + 7 (F) valence electrons, conserved
    // WHAT THIS GATES (re-derived 2026-07-23, post analytic-short/kappa/5-smooth + the 0.5(f2) raw-XC
    // feed): grid-continuation seeding makes the PRODUCTION fine grid converge CLEANLY to the aufbau
    // GROUND STATE.  Under raw-XC dynamics the XC-COLLAPSE basin is REMOVED (negCharge == 0 at every C in
    // the f2 acceptance sweep -- rho_DM >= 0 by construction), so this gate now carries basin history plus
    // the did-E-move pin.  The fine SCF converges in ~22 iters, charge conserved to 1e-8 throughout, at
    // -24.4325 -- 1.3 mHa from the CP2K SR2 truth -24.4312 (an Ecut=160-class number).  The historical
    // "-27.93 oracle / -3.5 Ha Exc step-2 gap" story recorded here previously was the RETRACTED
    // aliasing-era landscape (doc/GPWPlan.md TRAPS #2).
    EXPECT_TRUE(useSeed==false || scfF->Converged()) << "seeded fine SCF converges (no basin/spike thrash)";
    EXPECT_GT(Efine.GetTotalEnergy(), -29.0);   // basin avoidance: NOT the ~-40 unphysical attractor
    EXPECT_LT(Efine.GetTotalEnergy(), -20.0);   //                  NOT the +54 Pulay-thrash garbage
    EXPECT_NEAR(Efine.GetTotalEnergy(), -24.4304, 0.01);   // the raw-XC aufbau ground state at the production
                                                            //   default (auto Ecut=80, BallOnly); AliasFree@320
                                                            //   reference: -24.4325
}
