// File: IntegrationTests/GPW/Si.C  Si diamond -- the reference material (FCC a=10.26, 2-atom basis, Fd-3m non-symmorphic); the CP2K oracles live here.
//
// THE GRID (doc/TestSuitePlan.md §3): TEST(<Basis>_<Material>, <k>_[Grid]_[Fit]_[Sym]_[Spin]_[Occ]_[Machinery]_[Ansatz]_[Seed]_<Claim>)
// -- axis tokens in fixed order, the facade's defaults elided (k always named); the claim is CP2K (oracle anchor),
// Anchor (did-E-move), eq<Token> (a ONE-axis twin: this point equals the same point with that axis moved) or a
// property verb.  Tests are laid out in axis order.  `scripts/testgrid` renders the coverage table from these names.
//
//   GPW_Si.Γ_Imp_CP2K
//   GPW_Si.Γ_CP2K
//   GPW_Si.Γ_Imp_Schema
//   GPW_Si.Γ_RunsReal
//   GPW_Si.Γ_Deterministic
//   GPW_Si.Γ_eqCplx
//   GPW_Si.Γ_TranslationInvariant
//   GPW_Si.Γ_Imp_Smear_eqAufbau
//   GPW_Si.Γ_Imp_Pol_eqUnpol
//   GPW_Si.Γ_Imp_Pol_SpinSeed_eqUnpol
//   GPW_Si.Γ_Imp_M3_ShFermi_Smear_MomentRelaxes
//   GPW_Si.Γ_Imp_eqUnfolded
//   GPW_Si.Γ_Imp_Becke_eqUni
//   GPW_Si.Γ_Imp_Uni_DeltaFit_eqPWFit
//   GPW_Si.Γ_Imp_Uni_PWFit_Pol_eqUnpol
//   GPW_Si.k211_Imp_Anchor
//   GPW_Si.k222_CP2K
//   GPW_Si.k222_Imp_CP2K
//   GPW_Si.k222_Becke_Imp_eqFree
//   GPW_Si.k222s_Imp_CP2K
//   GPW_Si.k311_eqCplx
//   GPW_Si.k311_Uni_MOM_eqCplx

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


// (1) THE GAMMA ANCHOR == THE CP2K ENERGY GATE.  SR basis, Rcut=2a (every term translation-invariant and
// screened-complete), densityEcut=20 (FFT N=32): reproduces the CP2K FCC-Si Gamma GPW reference (SIPP_SR /
// GTH-PADE-q4 / LDA_X+VWN5) Etot=-7.11506 to the N=32 grid gap (~0.4 mHa; densityEcut>=30 -> -7.11505 exact).
// NOTE (analytic path): the old fast Rcut=0 anchor (-8.2476) is GONE -- the analytic collocation is always
// screened-complete Bloch, so home-only 1E matrices would MIX SCHEMES (Tr(D S_home)=8 while the grid density
// integrates the Bloch trace -- the forbidden inconsistency; doc/GPWPlan.md durable pins).  SR keeps the Bloch
// overlap cleanly PD at 2a.  Energy-gated at the density-fit floor (minDE=1e-6, minDrho relaxed to 1e-3).
TEST(GPW_Si, Γ_Imp_CP2K)
{
    const Material si=qchem::Materials::Get("Si_diamond");    // FCC primitive cell, 2-atom diamond basis, a=10.26
    const Lattice_3D lat=LatticeOf(si);
    SolidCalcOptions o=OptionsFor(si, "Si SR Gamma");
    o.densityEcut=20.0; o.imposeSymmetry=true;   // IMPOSED (the RunGPW default every Si anchor carried; V1.30)
    SCFParams par=ProductionGates();
    EnvOverrides(o, par);
    GpwReport report("Si "+o.label, par.Verbose);
    qchem::SolidCalculation calc(lat, MakeBasisSR(*si.cell), o, par);
    auto R=calc.Result();
    ASSERT_TRUE(R) << Why(R);
    EXPECT_NEAR(R->TotalCharge(), 8.0, 1e-6);      // 8 valence electrons
    EXPECT_NEAR(R->Energy(), -7.11506, 2e-3);      // CP2K FCC-Si Gamma reference (grid-gap tolerance)
}


//================================================================================================
TEST(GPW_Si, Γ_CP2K)
{
    const Material si=qchem::Materials::Get("Si_diamond");
    const Lattice_3D lat=LatticeOf(si);
    SolidCalcOptions o=OptionsFor(si, "Si SR Gamma (free, the facade's own recipe)");
    o.densityEcut=20.0;                              // FREE: the facade's default -- the imposed sibling is GPW_Si.Γ_Imp_CP2K
    qchem::SolidCalculation calc(lat, MakeBasisSR(*si.cell), o, ProductionGates());

    // N1/T1: the answers are reachable only through the PROOF, so a non-converged run cannot serve them.
    auto r = calc.Result();
    ASSERT_TRUE(r) << "SCF did not converge: " << (r ? std::string() : r.Error().details);
    EXPECT_NEAR(r->TotalCharge(), 8.0, 1e-6);
    EXPECT_NEAR(r->Energy(), -7.11506, 2e-3)        // the GPW_Si.Γ_Imp_CP2K anchor, same tolerance
        << "SolidCalculation must reproduce the driver's physics -- a facade that drifts from the "
           "recipe it fronts is worse than no facade";
    // The facade OWNS the grid decision, so it must be able to say what it chose (V1.26/V2.4): Si is soft
    // enough that the cost selector routes it to the uniform mesh, SIZED from this run's sharpness.
    EXPECT_EQ(calc.ResolvedXCMesh().cellKind, qcMesh::UnitCellKind::Uniform);
    EXPECT_GT(calc.ResolvedXCMesh().eCut, 0.0) << "an Auto-resolved uniform mesh must carry its own cutoff, "
                                                  "never fall back on nUniform's basis-blind default";
    // rho(r) is reachable through the same neutral ScalarFunction face the molecular facade exposes.
    EXPECT_GT(r->Density()(rvec3_t(0.1,0.1,0.1)), 0.0);
}


// The GPW run-report SCHEMA CHECK (RunReportPlan step 3).  Under an open run report the facade emits the
// `basis` (conditioning pre-flight) and `grids` (the ladder) sections itself during construction, and
// MakeIrrepWFs fills basis.perIrrep (per-Bloch-block conditioning) via the cursor.  Only the SETUP matters
// here, so a couple of iterations is plenty.
TEST(GPW_Si, Γ_Imp_Schema)
{
    const Material si=qchem::Materials::Get("Si_diamond");
    const Lattice_3D lat=LatticeOf(si);
    SolidCalcOptions o=OptionsFor(si, "Si grids schema");
    o.densityEcut=20.0; o.imposeSymmetry=true;

    report::ClearGlobal();                          // isolate this test's run
    {
        GpwReport report("Si "+o.label, false);
        qchem::SolidCalculation calc(lat, MakeBasisSR(*si.cell), o, Gates(4, 1e-3, 1e-6));
    }

    const report::json& all = report::GlobalReport();
    const report::json* grids=nullptr; const report::json* basis=nullptr;
    for (auto it=all.begin(); it!=all.end(); ++it)
    {
        if (it.value().contains("grids")) grids=&it.value()["grids"];
        if (it.value().contains("basis")) basis=&it.value()["basis"];
    }
    ASSERT_NE(grids, nullptr) << "no run emitted a grids section";
    EXPECT_GT((*grids)["densityEcut"].get<double>(), 0.0);
    EXPECT_TRUE(grids->contains("cutoffFactor"));
    EXPECT_TRUE(grids->contains("raster"));
    ASSERT_TRUE(grids->contains("ladder"));
    ASSERT_TRUE((*grids)["ladder"].is_array());
    ASSERT_GE((*grids)["ladder"].size(), 1u);
    const report::json& lvl0 = (*grids)["ladder"][0];
    for (const char* k : { "level", "N", "ecut", "nG", "role" })
        EXPECT_TRUE(lvl0.contains(k)) << "ladder row missing key: " << k;
    EXPECT_EQ(lvl0["role"].get<std::string>(), "reference");   // L==0 is the density/collocation grid
    ASSERT_TRUE(lvl0["N"].is_array());
    EXPECT_EQ(lvl0["N"].size(), 3u);
    ASSERT_TRUE(grids->contains("localPP"));
    EXPECT_TRUE((*grids)["localPP"].contains("kappa"));

    // The cursor path also gave GPW a basis section: per-Bloch-block conditioning in basis.perIrrep.
    ASSERT_NE(basis, nullptr) << "no run emitted a basis section";
    ASSERT_TRUE(basis->contains("perIrrep"));
    ASSERT_TRUE((*basis)["perIrrep"].is_array());
    ASSERT_GE((*basis)["perIrrep"].size(), 1u);
    EXPECT_TRUE((*basis)["perIrrep"][0].contains("cond"));
}


// The 3c-3 flip's REPORT EVIDENCE: basis.perIrrep now carries `runsReal` (the block's CONSTRUCTED
// scalar) beside `real` (the irrep's TRIM fact).  With the flip on every real:true block must actually
// RUN real -- the fact the acceptance twins above cannot see from outside; with it off (the harness
// default) the same blocks stay complex, so the two fields separate exactly where they should.
TEST(GPW_Si, Γ_RunsReal)
{
    auto perIrrep=[](){ const report::json& all=report::GlobalReport();
                        for (auto it=all.begin(); it!=all.end(); ++it)
                            if (it.value().contains("basis")) return it.value()["basis"]["perIrrep"];
                        return report::json{}; };
    const Material si=qchem::Materials::Get("Si_diamond");
    const Lattice_3D lat=LatticeOf(si);
    SolidCalcOptions o=OptionsFor(si, "Si real-block report");
    o.densityEcut=20.0;
    const SCFParams par=Gates(3, 1e-3, 1e-6);
    auto run=[&]{ GpwReport report("Si "+o.label, false); qchem::SolidCalculation calc(lat, MakeBasisSR(*si.cell), o, par); };

    report::ClearGlobal();
    o.forceComplex=false;                            // real TRIM blocks (the shipped default)
    run();
    report::json rows=perIrrep();
    ASSERT_GE(rows.size(), 1u);
    for (const auto& row : rows)
    {
        EXPECT_TRUE(row["real"].get<bool>());       // Γ-only: every block is TRIM
        EXPECT_TRUE(row["runsReal"].get<bool>());   // ...and the flip made each one ACTUALLY real
    }

    report::ClearGlobal();
    o.forceComplex=true;                             // the control: same irrep fact, complex build
    run();
    rows=perIrrep();
    ASSERT_GE(rows.size(), 1u);
    for (const auto& row : rows)
    {
        EXPECT_TRUE (row["real"].get<bool>());
        EXPECT_FALSE(row["runsReal"].get<bool>());
    }
}


//================================================================================================
//  CROSS-RUN DETERMINISM GATE (2026-08-18).  Three IDENTICAL SolidCalculation runs in one process
//  must give the SAME energy -- i.e. every run must replay fresh-process behaviour.  This is the
//  regression gate for the first-run anomaly: the GPW matrix-free 3C tensor closures capture the
//  evaluator's mutable per-SCF state (the CollocMemo D-screen), and when they were stored in the
//  process-wide DBCache a second identical run inherited the first run's converged D-screen -- its
//  seed Fock swept the full Hartree/Vxc pair set where a fresh process's (screened by the diagonal
//  SAD seed D) is diagonal-only: seed s-levels shifted ~0.24 Ha, converged E ~5e-6.  Fixed by
//  scoping those tensors per basis INSTANCE (tGPW_IBS::Overlap3C/Repulsion3C override); with the fix
//  the three energies here are BITWISE equal -- 1e-9 is pure headroom for BLAS/library variation.
//================================================================================================
TEST(GPW_Si, Γ_Deterministic)
{
    const Material si=qchem::Materials::Get("Si_diamond");
    const Lattice_3D lat=LatticeOf(si);
    const SCFParams par=ProductionGates();

    double E[3];
    for (int r=0;r<3;++r)
    {
        std::cout<<"[xrun] ================= RUN "<<r<<" ================="<<std::endl;
        qchem::SolidCalculation calc(lat, MakeBasisSR(*si.cell),
                                     {.Nelec=8, .species={{"Si",4}}, .densityEcut=20.0}, par);
        auto res=calc.Result();
        ASSERT_TRUE(res) << "run "<<r<<" did not converge: "<<res.Error().details;
        E[r]=res->Energy();
        std::cout.precision(15);
        std::cout<<"[xrun] run "<<r<<" E="<<E[r]<<std::endl;
    }
    std::cout.precision(15);
    std::cout<<"[xrun] E1-E0="<<E[1]-E[0]<<"  E2-E1="<<E[2]-E[1]<<std::endl;
    EXPECT_NEAR(E[1], E[0], 1e-9) << "cross-run pollution: first run differs";
    EXPECT_NEAR(E[2], E[1], 1e-9) << "runs 2+ should be steady";
}


//================================================================================================
//  STEP 3c-3 -- THE FACTORY-FLIP ACCEPTANCE (doc/RealComplexPlan.md).  SolidCalculation now computes
//  the working-type rule (irrep.IsReal() ∧ ham.PreservesReal()) and builds every TRIM block REAL by
//  default; forceComplex is the §6 ansatz-policy downgrade and this gate's A/B door.  The two runs
//  must be the SAME PHYSICS to machine precision: every per-term gate (RealComplexTerms.*) pinned the
//  real block's matrices bitwise except the quadrature GEMM's summation order (~1e-15/element), so
//  through a whole SCF the totals may differ only at accumulated-roundoff level.  A loose tolerance
//  here would be wrong twice over -- it could hide a genuine defect, and it would understate what the
//  per-term gates already guarantee.
//================================================================================================
//  TWO ARMS, because "machine-equal" means two different things across an SCF:
//  (1) ONE ITERATION from the shared uniform seed -- the PURE-ARITHMETIC gate.  Both runs build the
//      same Fock from the same D0 (statics + Hartree + raw XC bitwise per RealComplexTerms.*),
//      diagonalize, and fill the same gapped occupied set; every energy term from D1 must then agree
//      to roundoff.
//  (2) CONVERGED twins -- the PHYSICS gate.  Bitwise TRAJECTORIES are impossible by construction: the
//      real block diagonalizes with LAPACK's real path (dsyev) where the complex one runs zheev, which
//      agree only to roundoff even on an exactly-real matrix, and DIIS extrapolation chaotically
//      amplifies last-ulp seeds (measured: ~1e-4 in per-term energies after 12 iterations).  What the
//      physics guarantees is the FIXED POINT, so the converged states must agree to gate resolution.
//
//  (The 2026-08-18 WARM-UP DISCIPLINE is RETIRED: the first-run anomaly this gate used to hide from
//  -- a throwaway run so both arms sat in "steady-state slots" -- was the cross-run D-screen leak,
//  fixed by instance-scoping the GPW 3C tensors.  Every run now replays fresh-process behaviour
//  bit-for-bit, gated by GPW_Si.Γ_Deterministic above, so the arms run cold.)
static void ExpectRealComplexTwins(const Lattice_3D& lat, const Material& si, const char* what)
{
    const qchem::SolidCalcOptions optOn {.Nelec=8, .species={{"Si",4}}, .densityEcut=20.0};
    const qchem::SolidCalcOptions optOff{.Nelec=8, .species={{"Si",4}}, .densityEcut=20.0, .forceComplex=true};

    SCFParams one;                                   // the one-iteration recipe (gates 0 = strict-< never trips)
    one.NMaxIter=1; one.MinΔρ=0.0; one.MinΔE=0.0;
    one.MinΔFD=1e30; one.MinVirial=1e30; one.MinFD=1e30; one.StartingRelaxRo=0.3; one.MergeTol=1e-4;

    {   // (1) the arithmetic arm: exactly one iteration each
        qchem::SolidCalculation on(lat, MakeBasisSR(*si.cell), optOn, one), off(lat, MakeBasisSR(*si.cell), optOff, one);
        // ONE iteration by construction, so these are LAST-ITERATE diagnostics, not answers (N1/T1).
        const qchem::EnergyBreakdown En=on.LastIterateTerms(), Ec=off.LastIterateTerms();
        // 1e-8: real-vs-complex roundoff headroom -- Een's large-cancellation assembly (vs Enn~8 Ha)
        // measured 2e-9 across the 3-block mesh.  Roundoff grade, far below any defect scale.
        EXPECT_NEAR(En.GetTotalEnergy(), Ec.GetTotalEnergy(), 1e-8) << what << " (iteration 1)";
        EXPECT_NEAR(En["Kinetic"], Ec["Kinetic"], 1e-8) << what << " (iteration 1)";
        EXPECT_NEAR(En["Een"],     Ec["Een"],     1e-8) << what << " (iteration 1)";
        EXPECT_NEAR(En["Eee"],     Ec["Eee"],     1e-8) << what << " (iteration 1)";
        EXPECT_NEAR(En["Exc"],     Ec["Exc"],     1e-8) << what << " (iteration 1)";
        EXPECT_NEAR(on.LastIterateCharge(), off.LastIterateCharge(), 1e-10) << what << " (iteration 1)";
    }
    {   // (2) the physics arm: both converge on the production-shaped gates
        SCFParams par;
        par.NMaxIter=60; par.MinΔρ=1e-3; par.MinΔE=1e-6;
        par.MinΔFD=1e30; par.MinVirial=1e30; par.MinFD=1e30; par.StartingRelaxRo=0.3; par.MergeTol=1e-4;
        qchem::SolidCalculation on(lat, MakeBasisSR(*si.cell), optOn, par), off(lat, MakeBasisSR(*si.cell), optOff, par);
        auto ron=on.Result(), roff=off.Result();
        ASSERT_TRUE(ron)  << what << ": " << (ron  ? std::string() : ron .Error().details);
        ASSERT_TRUE(roff) << what << ": " << (roff ? std::string() : roff.Error().details);
        EXPECT_NEAR(ron->Energy(), roff->Energy(), 2e-5) << what;    // MinΔE=1e-6 resolution (measured 4.8e-6)
        const rvec3_t pts[]={ si.cell->ToCartesian(rvec3_t(0.3,0.4,0.7)),
                              si.cell->ToCartesian(rvec3_t(0.25,0.25,0.25)),
                              si.cell->ToCartesian(rvec3_t(0.1,0.9,0.2)) };
        for (const auto& r : pts)                                     // MinΔρ=1e-3 resolution (measured ≤9.1e-6)
            EXPECT_NEAR(ron->Density()(r), roff->Density()(r), 1e-4) << what;
    }
}

TEST(GPW_Si, Γ_eqCplx)
{
    const Material si=qchem::Materials::Get("Si_diamond");
    const Lattice_3D lat=LatticeOf(si);      // Γ-only: EVERY block is TRIM, so the flip makes the whole run real
    ExpectRealComplexTwins(lat, si, "Si Gamma real-vs-complex");
}


// DIAGNOSTIC: TERM-BY-TERM translation invariance of the 1E/PP matrix TRACES (no SCF -> fast) -- the tool
// that localized the Rcut>0 over-binding.  A rigid translation of the whole crystal (both atoms + their basis)
// must leave every trace invariant; the residual is that term's grid/mesh artifact.  Compare a CORNER atom
// (frac 0, on the cell boundary) vs an off-boundary atom (frac 0.13).  Kinetic is analytic -> the control.
//
// THE STORY (2026-07-09).  Pre-fix the KB nonlocal PP (Vnl) was translation-variant by ~16 Ha at ALL Rcut:
// MakeSeparablePP quadratured the RAW home orbital against the projector on a single-cell mesh, so a
// boundary-straddling corner orbital lost its wrapped tail (summing the PROJECTOR images cannot restore the
// ORBITAL's).  FIX = use the Bloch-summed orbital (Eval) as the bra (GPWPlan TODO 1a).  A NON-obvious extra:
// the local-PP/Hartree/XC (Vloc) variance was NOT the FFT raster (a uniform-grid ORIGIN shift is ~a no-op for
// a periodic quadrature -- Poisson summation only moves Nyquist-aliasing phases; the voxel-shift "Option A"
// was tried and REVERTED) -- it too was just incomplete ORBITAL WRAPPING.  Once the orbital is fully wrapped
// (Rcut>=2a) BOTH terms are translation-invariant to machine precision:
//     Rcut=1.50a  Kin d=0.0000  Vloc d=0.67    Vnl d=1.67
//     Rcut=2.00a  Kin d=0.0000  Vloc d=0.0000  Vnl d=0.0000   <-- fully wrapped
//     Rcut=3.00a  Kin d=0.0000  Vloc d=0.0000  Vnl d=0.0000
// (Committed Rcut=0 anchors are unaffected: at Rcut=0 itsRc={0} so Eval==the raw orbital.)  At Rcut=2a the
// residual is image-truncation limited (~1e-4, tightening with Rcut), so the guard tolerance is 1e-3 -- still
// 1600x below the ~1.7 Ha (pre-fix ~16 Ha) bug it protects against.
TEST(GPW_Si, Γ_TranslationInvariant)   // RE-ENABLED 2026-09-15: 0.3 s, invariant to 1e-10 -- it had no live reason to be parked
{
    using BasisSet::Gaussian::GPW_IBS;
    auto tr=[](const chmat_t& M){ double s=0; for (size_t i=0;i<M.rows();i++) s+=std::real(dcmplx(M(i,i))); return s; };
    const double a=10.26, dE=30.0;   // N=64 (finer than CP2K's converged grid)
    auto traces=[&](double frac, double& kin, double& vloc, double& vnl)
    {
        FCCUnitCell cell(a);
        cell.AddAtom(14, {frac, frac, frac});
        cell.AddAtom(14, {0.25+frac, 0.25+frac, 0.25+frac});
        Lattice_3D lat(cell, ivec3_t(1,1,1));
        auto st=lat.GetStructure();
        Pseudopotential::GTH_PP pp=Pseudopotential::GetGTH("Si","LDA",4);
        GPW_IBS gpw(cell, ivec3_t(1,1,1), ivec3_t(0,0,0), MakeBasisSR(cell), dE);  // eps-complete enumeration
        const BasisSet::Complex_OIBS& g=gpw;
        kin =tr(g.Kinetic());
        vloc=tr(gpw.MakeSpeciesFieldMatrix(st.get(), pp.local, qchem::BasisSet::FieldRange::Full));
        vnl =tr(gpw.MakeProjectorMatrix(st.get(), pp.nonlocal));
    };
    {
        double kc,lc,nc, ks,ls,ns;
        traces(0.00, kc,lc,nc);
        traces(0.13, ks,ls,ns);
        std::printf("Kin[%10.5f/%10.5f d=%.2e]  Vloc[%10.5f/%10.5f d=%.2e]  Vnl[%10.5f/%10.5f d=%.2e]\n",
                    kc,ks,std::fabs(kc-ks), lc,ls,std::fabs(lc-ls), nc,ns,std::fabs(nc-ns));
        EXPECT_NEAR(lc, ls, 1e-3) << "Vloc translation invariance (complete enumeration)";
        EXPECT_NEAR(nc, ns, 1e-3) << "Vnl translation invariance (complete enumeration)";
    }
}


// (4b-i) FERMI SMEARING IS INERT ON A GAP (doc/GPWPlan1.md 4b, gate i).  The same gapped Si/Gamma anchor as
// GPW_Si.Γ_Imp_CP2K, but with smearing kT=1e-3 Ha turned ON.  Si is a wide-gap insulator in this basis
// (ε_LUMO−ε_HOMO ≫ kT), so every f_i is essentially 0 or 1: the fractional occupations collapse to the
// aufbau integers, the Mermin −TS is negligible, and the total reproduces the CP2K reference −7.11506.  This
// is the regression that smearing must not perturb a system that does not need it (the T→0 / kT≪gap limit).
TEST(GPW_Si, Γ_Imp_Smear_eqAufbau)
{
    const Material si=qchem::Materials::Get("Si_diamond");
    const Lattice_3D lat=LatticeOf(si);
    SolidCalcOptions o=OptionsFor(si, "Si SR Gamma +smear");
    o.densityEcut=20.0; o.imposeSymmetry=true;
    o.seed=qchem::ChargeDensity::SeedStrategy::Uniform;
    SCFParams par=ProductionGates();
    par.SmearingkT=1e-3;
    EnvOverrides(o, par);
    GpwReport report("Si "+o.label, par.Verbose);
    qchem::SolidCalculation calc(lat, MakeBasisSR(*si.cell), o, par);
    auto R=calc.Result();
    ASSERT_TRUE(R) << Why(R);
    EXPECT_NEAR(R->TotalCharge(), 8.0, 1e-6);
    EXPECT_NEAR(R->Energy(), -7.11506, 2e-3);                // == the no-smear anchor: smearing is inert on a gap
    EXPECT_NEAR(R->EnergyTerms()["MinusTS"], 0.0, 1e-4);     // −TS negligible when kT ≪ gap (f_i ∈ {0,1})
}


// (tier 4b, invariant) THE ζ=0 COLLAPSE: the TWO-CHANNEL machinery on a CLOSED SHELL must reproduce the
// unpolarized anchor.  Same gapped Si/Γamma cell + recipe as GPW_Si.Γ_Imp_Smear_eqAufbau, but multiplicity=1 drives
// the polarized pipeline (the dcmplx composite WF under SpinGroup::Polarized, Crystal_EC(4,4), the spin-native XC term) with
// nUp=nDn=4 -- v^σ(ρ/2,ρ/2)=v^P(ρ) pointwise, so the total must land on the SAME −7.11506 anchor.  The
// periodic sibling of the molecular WaterPolarizedLDA-vs-LDA check; catches any polarized-path divergence
// (channel bookkeeping, shared-engine caching, the collocation memo screen) on known ground.
TEST(GPW_Si, Γ_Imp_Pol_eqUnpol)
{
    const Material si=qchem::Materials::Get("Si_diamond");
    const Lattice_3D lat=LatticeOf(si);
    SolidCalcOptions o=OptionsFor(si, "Si SR Gamma pol-singlet");
    o.imposeSymmetry=true;   // V1.30: was the DEFAULT; now stated, because an imposition you did not ask for is invisible in the result
    o.multiplicity=1;                                  // EXPLICIT two-channel singlet (nUp=nDn=4)
    o.densityEcut=20.0;
    GpwReport report("Si "+o.label, false);
    qchem::SolidCalculation calc(lat, MakeBasisSR(*si.cell), o, ProductionGates());
    auto R=calc.Result();
    ASSERT_TRUE(R) << Why(R);
    EXPECT_NEAR(R->TotalCharge(), 8.0, 1e-6);
    EXPECT_NEAR(R->Energy(), -7.11506, 2e-3);          // == the unpolarized Becke anchor (ζ=0 collapse exact)
}


TEST(GPW_Si, Γ_Imp_Pol_SpinSeed_eqUnpol)
{
    const Material si=qchem::Materials::Get("Si_diamond");
    const Lattice_3D lat=LatticeOf(si);
    SolidCalcOptions o=OptionsFor(si, "Si SR Gamma pol-singlet spin-SAD");
    o.imposeSymmetry=true;   // V1.30: was the DEFAULT; now stated, because an imposition you did not ask for is invisible in the result
    o.multiplicity=1;                                  // EXPLICIT two-channel singlet (nUp=nDn=4)
    o.densityEcut=20.0;
    o.seed=qchem::ChargeDensity::SeedStrategy::SAD;    // -> PolarizedSeedCD (rho/2 channels for pairless Si)
    GpwReport report("Si "+o.label, false);
    qchem::SolidCalculation calc(lat, MakeBasisSR(*si.cell), o, ProductionGates());
    auto R=calc.Result();
    ASSERT_TRUE(R) << Why(R);
    EXPECT_NEAR(R->TotalCharge(), 8.0, 1e-6);
    EXPECT_NEAR(R->Energy(), -7.11506, 2e-3);          // the SAME unpolarized Becke anchor
}


// The SPIN-SAD sibling of the ζ=0 collapse (SCFSeedingPlan §10 increment B): a polarized run with a SAD
// seed now assembles the TWO-CHANNEL PolarizedSeedCD -- Si's library entry is spin-agnostic (closed shell),
// so each channel is exactly rho/2 and the whole polarized-seed machinery (channel SeedCDs, the merged
// FourierDensity total into Hartree, the cSpinResolved_CD branch in RhoPol) must land on the SAME anchor.
// The proof the polarized seed cannot perturb non-magnetic physics.
// ONE CHEMICAL POTENTIAL OVER BOTH SPIN CHANNELS (Crystal_EC spinsShareFermi; 2026-08-10).
//
// WHAT IT PINS: with a shared reservoir the MOMENT IS AN OUTPUT.  Seed Si-Gamma -- a closed-shell,
// gapped, non-magnetic system -- as a TRIPLET (nUp=5, nDn=3) and let the fill decide: one mu must move the
// two excess-spin electrons back down and land on the SINGLET energy.  With separate per-channel counts
// (the control below) nUp-nDn=2 is conserved by construction, so the run is stuck on the triplet no matter
// how far it converges.  The difference between the two arms IS the ensemble.
//
// WHY IT MATTERS (doc/SymmetryUpgradePlan.md sec 7 step 7): separate reservoirs mean mu_up != mu_dn, and then
// an occupation is monotone in epsilon only WITHIN a channel -- MnO run 29 ended with a down level 27 mHa
// BELOW an up level and LESS occupied.  That gap is an unrelieved driving force to move charge between the
// channels, so the converged state is not the free minimum.  CP2K's MnO deck constrains nothing.
TEST(GPW_Si, Γ_Imp_M3_ShFermi_Smear_MomentRelaxes)
{
    const Material si=qchem::Materials::Get("Si_diamond");
    const Lattice_3D lat=LatticeOf(si);
    struct Arm { double charge, E; };
    auto TripletSi=[&](bool shared) -> Arm
    {
        SolidCalcOptions o=OptionsFor(si, shared ? "Si Gamma triplet, shared mu" : "Si Gamma triplet, per-channel mu");
        o.imposeSymmetry=true;
        o.multiplicity=3;                              // SEEDED as a triplet (nUp=5, nDn=3)
        o.densityEcut=20.0;
        o.spinsShareFermi=shared;
        SCFParams par=ProductionGates();
        par.SmearingkT=5e-3;                           // a shared mu needs smearing to relax the moment with
        GpwReport report("Si "+o.label, false);
        qchem::SolidCalculation calc(lat, MakeBasisSR(*si.cell), o, par);
        // Both arms are read as LAST ITERATES on purpose: the held arm is a constrained state that need not
        // meet the production gates, and the comparison below is between the two arms, not against a pin.
        return {calc.LastIterateCharge(), calc.LastIterateTerms().GetTotalEnergy()};
    };

    const Arm shared=TripletSi(true);
    EXPECT_NEAR(shared.charge, 8.0, 1e-6);             // the TOTAL is what a shared reservoir conserves
    // The moment relaxed away: this lands on the singlet anchor of GPW_Si.Γ_Imp_Pol_eqUnpol.
    EXPECT_NEAR(shared.E, -7.11506, 3e-3)
        << "a shared mu must let a triplet-seeded closed-shell system fall back to the singlet";

    // CONTROL: the same run with separate per-channel counts cannot relax -- nUp-nDn=2 is conserved.
    const Arm held=TripletSi(false);
    EXPECT_NEAR(held.charge, 8.0, 1e-6);
    EXPECT_GT(held.E, shared.E + 1e-3)
        << "with two reservoirs the seeded multiplicity is a CONSTRAINT and must sit above the free minimum";
}


// (T3.2, doc/SymmetryUpgradePlan.md §6b) STREAM FOLD through-SCF A/B on an IMPOSED Γ-only run: the factory
// arms route (b) on the shared molecular evaluator (reduced stream build + replay, rep-transform h), and the
// existing SymmetrizeGMap/SymmetrizeRaster sites complete the group-average.  GPW_STREAM_FOLD=0/1 toggles the
// fold between two otherwise identical runs (read fresh in the factory, so one process can A/B).  The two
// totals must agree to the band-limit class (§8 through-SCF tier -- the production 5-smooth grid is NOT
// τ-commensurate, so reduced+P and full+P are two equally valid quadratures of the same density), and the
// folded run's [stream cache] line must show the reduced build (repPairs << pairs).
TEST(GPW_Si, Γ_Imp_eqUnfolded)
{
    const Material si=qchem::Materials::Get("Si_diamond");    // Fd-3m: non-symmorphic, 48 ops, quarter glide
    const Lattice_3D lat=LatticeOf(si);                       // Γ-only: the T3.2 arming condition (k≠Γ is T3.4)
    SolidCalcOptions o=OptionsFor(si, "Si diamond Γ stream-fold A/B");
    o.densityEcut=20.0; o.imposeSymmetry=true;
    o.seed=qchem::ChargeDensity::SeedStrategy::Uniform;
    struct Arm { double charge, E; };
    auto arm=[&](const char* fold) -> Arm
    {
        setenv("GPW_STREAM_FOLD",fold,1);  qchem::ReresolveRunPolicy();
        GpwReport report("Si "+o.label, false);
        qchem::SolidCalculation calc(lat, MakeBasisSR(*si.cell), o, ProductionGates());
        auto R=calc.Result();
        EXPECT_TRUE(R) << Why(R);
        return {calc.LastIterateCharge(), calc.LastIterateTerms().GetTotalEnergy()};
    };
    const Arm R0=arm("0"), R1=arm("1");
    unsetenv("GPW_STREAM_FOLD");  qchem::ReresolveRunPolicy();

    std::cout << "[stream fold A/B] E(full)=" << R0.E << "  E(folded)=" << R1.E << "  dE=" << R1.E-R0.E << std::endl;
    EXPECT_NEAR(R1.E, R0.E, 1e-5)
        << "route (b) reduced streams must reproduce the full-stream imposed Γ run (band-limit class)";
    EXPECT_NEAR(R1.charge, 8.0, 1e-6);
}


//================================================================================================
//  THE BECKE XC GATE (doc/GPWPlan1.md "Becke XC grid").  The atom-centred periodic Becke XC
//  quadrature must reproduce the uniform-multigrid XC on a CONDITIONED basis -- same E_xc and same
//  V_xc matrix to grid tolerance -- before it can become the default for diffuse bases.  Converge
//  Si/Gamma on the standard uniform route (the GPW_Si.Γ_Imp_CP2K recipe), then evaluate BOTH
//  XC term pairs (Dirac + VWN5) on the SAME converged density:
//    uniform -- the PAIR quadrature on the Vxc fit basis's FFT grid (the raw-collocation route);
//    Becke   -- Vxc_Quadrature: rho(r) analytic at the atom-centred points, MatrixOverlap matrix.
//  Angular rule: GaussLegendre (machine-exact algebraic degree at any L -- the audited Lebedev
//  tables stop at L=11; see the Mesh_AngularDegree tests).
//================================================================================================
TEST(GPW_Si, Γ_Imp_Becke_eqUni)
{
    const Material si=qchem::Materials::Get("Si_diamond");
    const Lattice_3D lat=LatticeOf(si);

    // Converge on the standard uniform route.  densityEcut=60 (not the SCF-sufficient 20): the gate
    // compares MATRIX ELEMENTS, and the uniform raw-adjoint H_xc carries raster error at N=15^3 that the
    // ~1e-4-converged Becke quadrature exposes (measured at Ecut=20: dExc=1.2e-4 but max|U-B|=1.2e-2 --
    // the raster's error, not Becke's; at Ecut=60 max|U-B|=3.5e-4).
    SolidCalcOptions o=OptionsFor(si, "Si Becke gate");
    o.imposeSymmetry=true;   // V1.30: was the DEFAULT; now stated, because an imposition you did not ask for is invisible in the result
    o.densityEcut=60.0;
    o.xcMesh.cellKind=qcMesh::UnitCellKind::Uniform;   // the gate's SCF arm is the uniform route BY DESIGN (Auto would flip it)
    SCFParams par=ProductionGates(); par.MergeTol=SCFParams{}.MergeTol;   // (this gate never set MergeTol)
    GpwReport report("Si "+o.label, false);
    qchem::SolidCalculation calc(lat, MakeBasisSR(*si.cell), o, par);
    ASSERT_TRUE(calc.Result()) << Why(calc.Result());
    const GpwHandles h=Handles(calc);

    XCProbe U=UniformXCProbe(h, lat.GetStructure());
    XCProbe B=BeckeXCProbe(h, lat.GetStructure(), "B40", qcMesh::BeckeXCParams());
    EXPECT_NEAR(B.Exc, U.Exc, 5e-4);                 // measured: dExc=1.1e-4
    EXPECT_NEAR(B.rhoLost, 0.0, 5e-3);               // the Becke mesh integrates rho to Tr(DS)
    EXPECT_LT(DiffXC(U,B), 1e-3);                    // measured: 3.5e-4
}


// ===== EXPERIMENTAL (scratch): global-μ across k-blocks (item 3 inc 3) =====
// AL_KGRID=n (mesh nxnxn), AL_GLOBAL=0/1 (per-block vs global μ), AL_KT, AL_NMAX.
// I2 (plan §6a fit/grid SEPARATION): the (Delta, uniform) cross cell.  The SAME material through the
// PLANE-WAVE fit (band-limited v_xc on the FFT raster) and through the DELTA fit on the uniform cell
// mesh -- the two v_xc representations must agree to the route-difference class (band-limiting +
// raster-vs-midpoint-mesh quadrature), the same class the Becke-vs-uniform gate measures (~1e-4 Exc).
TEST(GPW_Si, Γ_Imp_Uni_DeltaFit_eqPWFit)
{
    const Material si=qchem::Materials::Get("Si_diamond");
    const Lattice_3D lat=LatticeOf(si);
    SolidCalcOptions o=OptionsFor(si, "Si PW-fit");
    o.imposeSymmetry=true;   // V1.30: was the DEFAULT; now stated, because an imposition you did not ask for is invisible in the result
    o.densityEcut=20.0;
    o.seed=qchem::ChargeDensity::SeedStrategy::Uniform;
    SCFParams par=ProductionGates(); par.MergeTol=SCFParams{}.MergeTol;
    auto energy=[&]{ GpwReport report("Si "+o.label, false);
                     qchem::SolidCalculation calc(lat, MakeBasisSR(*si.cell), o, par);
                     auto r=calc.Result(); EXPECT_TRUE(r) << Why(r);
                     return r ? r->Energy() : 0.0; };

    o.xcMesh.cellKind=qcMesh::UnitCellKind::Uniform;                        // (PlaneWave, raster)
    const double P=energy();

    o.label="Si Delta-fit uniform"; o.vxcFit=Hamiltonian::VxcFit::Delta;    // (Delta, uniform mesh)
    o.xcMesh.eCut=o.densityEcut;                                            // resolve rho on the midpoint mesh
    const double D=energy();

    EXPECT_NEAR(D, P, 5e-3)
        << "the delta fit on the uniform cell mesh must reproduce the PW fit on the raster to the "
           "band-limit/quadrature route-difference class (plan 6a fit/grid separation)";
}


// THE V2.3 GATE (doc/CleanupCandidates.md V2.3, run 2026-09-14): the POLARIZED PLANE-WAVE fit route.  The row
// said a polarized Ham_PW_DFT THROWS under VxcFit::PlaneWave because the raster route had no per-spin rho
// cache; the raster route grew RhoPol/RefreshPol on 2026-08-28 and since V1.37 step 3 the one XC term asks it
// for the pair on any fit basis -- but nothing had ever RUN it.  Same cell and recipe as the (PlaneWave,
// raster) arm above, as the EXPLICIT two-channel singlet (nUp=nDn=4): the zeta=0 collapse must land on the
// unpolarized PW-fit answer, exactly as GPW_Si.Γ_Imp_Pol_eqUnpol pins it on the Becke route.
TEST(GPW_Si, Γ_Imp_Uni_PWFit_Pol_eqUnpol)
{
    const Material si=qchem::Materials::Get("Si_diamond");
    const Lattice_3D lat=LatticeOf(si);
    SolidCalcOptions o=OptionsFor(si, "Si PW-fit unpol");
    o.imposeSymmetry=true;
    o.densityEcut=20.0;
    o.seed=qchem::ChargeDensity::SeedStrategy::Uniform;
    o.xcMesh.cellKind=qcMesh::UnitCellKind::Uniform;   // (PlaneWave, raster): Auto resolves the PW fit here
    SCFParams par=ProductionGates(); par.MergeTol=SCFParams{}.MergeTol;
    struct Arm { double charge, E; };
    auto arm=[&]() -> Arm { GpwReport report("Si "+o.label, false);
                            qchem::SolidCalculation calc(lat, MakeBasisSR(*si.cell), o, par);
                            auto r=calc.Result(); EXPECT_TRUE(r) << Why(r);
                            return r ? Arm{r->TotalCharge(), r->Energy()} : Arm{0,0}; };
    const Arm U=arm();

    o.label="Si PW-fit pol-singlet"; o.multiplicity=1;   // the explicit two-channel singlet on the SAME route
    const Arm P=arm();
    EXPECT_NEAR(P.charge, 8.0, 1e-6);
    EXPECT_NEAR(P.E, U.E, 1e-6)                            // measured 6e-9 on 2026-09-14, same 17 iterations
        << "zeta=0 collapse on the plane-wave fit route: v^sigma(rho/2,rho/2) == v(rho) pointwise, so the "
           "two-channel singlet must reproduce the unpolarized PW-fit answer";
}


// (1) THE REAL-MATERIAL SCF: crystalline silicon (diamond) primitive cell at Gamma, driven end-to-end by
// the framework cSCFIterator through the plane-wave Kohn-Sham Hamiltonian on a GAUSSIAN (GPW) basis.  8
// valence electrons (2 x Zion 4) fill a closed shell (sigma_g^2 sigma_u^2 pi_u^4), so it converges cleanly.
// With the G-space local PP (box-independent, PW G=0/alignment convention) the total is now PHYSICAL:
// Etot=-8.248 -- close to the plane-wave bulk-Si -7.2273 (Ecut=4) / converged ~-7.9.  The residual ~1 Ha is
// the Rcut=0 over-binding (home-cell electrons, no inter-cell screening, feel the full periodic ion Ewald);
// true bulk (Rcut>0) awaits the overlap-conditioning fix.  A did-E-move regression anchor (pin the value).
// (1c) MULTI-K PLUMBING: a 2x1x1 Monkhorst-Pack mesh (2 k-points), SR basis, Rcut=2a.  The ANALYTIC
// collocation always sums the screened cross-cell pair offsets (with their Bloch phases), so there is no
// "Rcut=0 makes every k-block a copy of Gamma" shortcut any more -- the 2-point mesh has REAL dispersion, and
// this pins its BZ-weighted total.  What the gate protects is the multi-k machinery: one GPW_IBS per BZ
// k-point (GPW_BasisSet iterating MakeKMesh WITH BZ weights), the multi-block GetIrreps, Crystal_EC's
// BZ-weighted (Sum_k w_k) occupation, the per-irrep k-loop, and the BZ-summed charge/energy (it caught a
// missing BZ weight -> charge x Nk).  Energy-gated at the fit floor like the Gamma anchor.
// ★ RE-ENABLED 2026-09-09, AND IT IS THE GATE ON KP-0.  Disabled since 2026-07-20; when it was run again
// on 2026-09-09 it reported charge=12 against 8 valence electrons, with its own banner saying why --
// "[IBZ] 2 k-points -> 2 irreducible; Σw=1.5", and 8 x 1.5 = 12 exactly.  The IBZ star weights did not
// sum to 1 because `FoldGrid` applied the reciprocal op U to the grid INDEX vector and let the mod-N wrap
// place the image: on an ANISOTROPIC 2x1x1 mesh an axis-permuting cubic op is not a mesh symmetry at all,
// so the action was not a bijection and the two "orbits" overlapped (sizes 1 and 2 on a 2-point grid).
// The index action is the CONJUGATED matrix M=D U D^{-1}; ops with non-integral M are not mesh symmetries
// and are now skipped wholesale (src/Symmetry/Lattice_3D/Imp/Fold.C, gates in src/Symmetry/tests/L_Fold.C).
// The charge assertion below is the one that caught it, and it is the kind that CANNOT go stale: 8 valence
// electrons is physics, not a banked number.  It runs in ~3 s.
TEST(GPW_Si, k211_Imp_Anchor)
{
    const Material si=qchem::Materials::Get("Si_diamond");
    const Lattice_3D lat=LatticeOf(si, ivec3_t(2,1,1));   // 2 k-points: Gamma + the zone-boundary k=1/2 (real +-1 phases)
    SolidCalcOptions o=OptionsFor(si, "Si SR 2x1x1");
    o.densityEcut=20.0; o.imposeSymmetry=true;   // IMPOSED, as every RunGPW anchor was (its default; V1.30)
    SCFParams par=ProductionGates();
    EnvOverrides(o, par);
    GpwReport report("Si "+o.label, par.Verbose);
    qchem::SolidCalculation calc(lat, MakeBasisSR(*si.cell), o, par);
    auto R=calc.Result();
    ASSERT_TRUE(R) << Why(R);
    EXPECT_NEAR(R->TotalCharge(), 8.0, 1e-6);                  // 8 valence e- (BZ-weighted Sum_k, not x Nk)
    // Did-E-move anchor, RE-PINNED 2026-09-09 to the value measured with correct weights (was -7.45137,
    // banked 2026-07-15 before the IBZ fold path existed; the run drifted 1.6e-3 over the intervening GPW
    // work and stayed inside this 5e-3 fit-floor window throughout).  JUDGED, not just refreshed: the
    // Gamma-point 2x1x1 SUPERCELL -- band-folding-equivalent, a fully independent route -- gives
    // -7.451621 per primitive cell (DISABLED_SiSupercellLadder, SI_LADDER=2,1,1, itself still descending
    // at its iteration cap), so the two routes agree to 1.3e-3 and bracket the old anchor.
    EXPECT_NEAR(R->Energy(), -7.45294, 5e-3);
}


// DISPERSIVE MULTI-K BULK -- RE-ENABLED 2026-09-15 (5.8 s after the 2026-08 box-walk work; it was parked at
// "~4 min" in July) -- the first REAL bulk GPW, unblocked by the KB
// Bloch-orbital fix (Rcut>0 now correct).  Gamma-centred 2x2x2 MP, SIPP_SR, Rcut=2a: charge stays 8 and the
// total drops with k-sampling (Gamma -7.11467 -> 2x1x1 -7.451 -> 2x2x2 -7.778 -- real dispersion).
// CROSS-CHECK vs CP2K AT THE SAME GAMMA-CENTRED MESH: -7.7778 vs CP2K -7.77846 (~0.7 mHa, the N=32 grid gap;
// deck UnitTests/CP2K/si_fcc_gpw_222_gamma.inp).  The 90 mHa vs CP2K's DEFAULT -7.86744 is the k-CONVENTION:
// Gamma-centred here (kShift=0) vs CP2K's classic SHIFTED MONKHORST-PACK (k at +/-1/4).  The shifted grid is
// the sibling test below (kShift=1/2).  The general-k PHYSICS is validated at both.
TEST(GPW_Si, k222_CP2K)
{
    const Material si=qchem::Materials::Get("Si_diamond");
    const Lattice_3D lat=LatticeOf(si, ivec3_t(2,2,2));
    SolidCalcOptions o=OptionsFor(si, "Si 2x2x2 Gamma-centred (free)");
    o.densityEcut=20.0; o.imposeSymmetry=false;      // FREE (2026-09-15): the full 8-k mesh; GPW_Si.k222_Imp_CP2K is the IMPOSED arm at the same anchor
    SCFParams par=TightGates(60);
    EnvOverrides(o, par);
    GpwReport report("Si "+o.label, par.Verbose);
    qchem::SolidCalculation calc(lat, MakeBasisSR(*si.cell), o, par);
    auto R=calc.Result();
    ASSERT_TRUE(R) << Why(R);
    EXPECT_NEAR(R->TotalCharge(), 8.0, 1e-6);
    EXPECT_NEAR(R->Energy(), -7.77846, 3e-3) << "GPW 2x2x2 Gamma-centred vs CP2K same-mesh -7.77846";
}


// (item 5, IBZ) NON-SYMMORPHIC -- diamond Si (FCC lattice + a 2-atom basis at (0,0,0),(¼,¼,¼)) is space group
// Fd-3m: NON-symmorphic (the two sublattices are related by a glide, τ=(¼,¼,¼)≠0).  The k-FOLD reaches the full
// Oh (3 irreducible k-points, same as FCC Al on this lattice): the τ=0 Td subgroup + time reversal (k→−k) already
// supplies the inversion Td lacks.  The DENSITY star-average now carries the glide τ: the G-space Hartree via the
// e^{+2πi(Um)·τ} phase (SymmetrizeGMap over SpaceGroup::ReciprocalOps) and the real-space XC raster via the exact
// FFT fractional shift ρ(W·x+τ) (SymmetrizeRaster over SpaceGroup::DirectOps).  So the IBZ-reduced total now
// reproduces the full-mesh Γ-centred 2×2×2 exactly (was −8.259, ~0.48 Ha off, under the old τ=0 W-only guard).
TEST(GPW_Si, k222_Imp_CP2K)
{
    const Material si=qchem::Materials::Get("Si_diamond");    // diamond = FCC + the glide-related 2nd sublattice (non-symmorphic)
    const Lattice_3D lat=LatticeOf(si, ivec3_t(2,2,2));
    SolidCalcOptions o=OptionsFor(si, "Si diamond IBZ");
    o.densityEcut=20.0; o.imposeSymmetry=true;
    o.seed=qchem::ChargeDensity::SeedStrategy::Uniform;
    GpwReport report("Si "+o.label, false);
    qchem::SolidCalculation calc(lat, MakeBasisSR(*si.cell), o, ProductionGates());
    auto R=calc.Result();
    ASSERT_TRUE(R) << Why(R);
    // Target: the IBZ-reduced total reproduces the full-mesh Γ-centred 2×2×2 (DISABLED_SR_2x2x2GammaCentred,
    // -7.77846) to grid/SCF tolerance -- the non-symmorphic glide τ-phase makes the reduced density exact
    // (measured -7.77847, ~1e-5 vs the full mesh; the fold reaches 3 irreducible k-points under the full Oh).
    EXPECT_NEAR(R->Energy(), -7.77846, 2e-3)
        << "diamond Si (non-symmorphic Fd-3m): IBZ density symmetrization with the glide τ-phase must match the "
           "full mesh -7.77846 (G-space e^{+2πi(Um)·τ} + real-space FFT τ-shift)";
}


// THE W1 GATE (doc/SymmetryUpgradePlan.md §6a): Becke XC under IBZ.  BeckeFit_IBS group-averages its
// mesh INVARIANT and star-averages rho every iteration (the SymmetrizeRaster hook, exact orbit-mean
// projector); on the converged SYMMETRIC density the invariant mesh integrates identically to the
// single-orientation mesh (Q_inv(f)==Q(f) for a symmetric f), so the IBZ run must reproduce the
// full-mesh Becke run to the IBZ class -- the Becke sibling of GPW_Si.k222_Imp_CP2K, with the
// non-symmorphic glide tau exercised through the torus fold + MakeInvariant.  COARSE explicit Becke
// recipe on BOTH arms: the gate compares like against like, so grid quality cancels (and the
// imposed arm's group-average mesh growth stays affordable).
TEST(GPW_Si, k222_Becke_Imp_eqFree)
{
    const Material si=qchem::Materials::Get("Si_diamond");
    const Lattice_3D lat=LatticeOf(si, ivec3_t(2,2,2));
    SolidCalcOptions o=OptionsFor(si, "Si diamond Becke FULL");
    o.densityEcut=20.0;
    o.xcMesh=qcMesh::BeckeXCParams(15, 2.0, 9);          // explicit coarse Becke (nR=15, GL-9), same on both arms
    o.seed=qchem::ChargeDensity::SeedStrategy::Uniform;
    auto energy=[&]{ GpwReport report("Si "+o.label, false);
                     qchem::SolidCalculation calc(lat, MakeBasisSR(*si.cell), o, ProductionGates());
                     auto r=calc.Result(); EXPECT_TRUE(r) << Why(r);
                     return r ? r->Energy() : 0.0; };

    o.imposeSymmetry=false;
    const double F=energy();
    o.label="Si diamond Becke IBZ"; o.imposeSymmetry=true;
    const double R=energy();

    std::cout<<"[Becke IBZ gate] full="<<F<<" reduced="<<R<<" dE="<<R-F<<std::endl;
    // Tolerance = the RULE-DIFFERENCE class at this deliberately coarse L: the free arm runs GL-9
    // (50 dirs), the imposed arm the MIXED-rule site-adapted minimal grid (degree-9-exact, ~76
    // dirs/atom) -- measured 2.0e-3 at L=9, collapsing with L (passes 2e-3 already at L=17; both
    // rules sit on the comparison floor at the production L=29).
    EXPECT_NEAR(R, F, 3e-3)
        << "Becke+IBZ must reproduce Becke+full-mesh (the W1 star-average makes the reduced density exact "
           "on the invariant Becke mesh; doc/SymmetryUpgradePlan.md 6a)";
}


// SHIFTED Monkhorst-Pack (kShift=½ → k at ±¼ = CP2K's DEFAULT MONKHORST-PACK 2 2 2) -- the apples-to-apples
// match to CP2K's shipped 2x2x2 reference -7.86744 (deck si_fcc_gpw_222.inp).  This is the FIRST run with a
// genuinely COMPLEX Bloch phase e^{ik·R} (not ±1), so the density matrix D and every k-block matrix are
// genuinely complex -- the exact case that exposed (and now validates the fix for) TWO complex-only bugs
// (doc/GPWPlan.md, "Complex-k GPW FIXED", 2026-07-10):
//   1. GPW_Evaluator::BuildWeights conjugated the BRA (i) instead of the KET (j) collocation slot, so the
//      Fourier density rho-tilde was the TRANSPOSE-density D^T (a different real field at complex k) -> the
//      Hartree/XC drive was inconsistent with the physical density (IrrepCD::operator() / the PW delta path).
//   2. GPW_Evaluator::MakeSeparablePP summed the KB projector images with e^{+ik·R} instead of e^{-ik·R};
//      the correct Bloch projection b_i=<chi_i^k|beta_home> tiles all-space with a CONJUGATED image phase.
//      At complex k this HALVED the nonlocal-PP trace (Vnl 42->22) -> a spurious deep core level -> over-bind.
// Both are inert at Gamma / half-integer k (phase ±1 self-conjugate), so every committed anchor is unchanged;
// they matter ONLY here.  Grid-matched to CP2K's mesh.
// ENABLED 2026-07-15 (\S0a complex-k revalidation): the ANALYTIC collocate/integrate kernels reproduce the
// CP2K shifted-MP reference -- Rcut=2a gave -7.86724 (0.20 mHa), and the run is affordable now (the stream
// cache + the phase-independent integrate memo make the 8 k-blocks share the static sweeps: ~2.5 min).
// Rcut switched to AUTO for scheme consistency with the enabled anchors (both sides parameter-free).
//
// RE-ENABLED (undisabled) 2026-08-19, and it is the point of the whole exercise: while this test sat
// DISABLED it silently ROTTED to -3.7351, and nothing caught it because it is the ONLY fractional-k SCF
// coverage in the suite -- every other k in every other test is TRIM, where the Bloch phases are +-1 and
// the defect is structurally invisible.  The cause was the D-aware integrate-back screen testing
// |Re(D_ij e^{-ik.R_n})| as if it were a magnitude: at a quarter-integer k that real part vanishes for
// every ODD offset, so the Hartree/XC matrix lost its entire imaginary part (doc/Benchmark.md footnote 1).
// It now runs in ~14 s (not the ~2.5 min the note above records), so there is no cost argument for hiding
// it again.  A DISABLED regression test is a test that will be wrong when you next need it.
TEST(GPW_Si, k222s_Imp_CP2K)
{
    const Material si=qchem::Materials::Get("Si_diamond");
    const Lattice_3D lat=LatticeOf(si, ivec3_t(2,2,2));
    SolidCalcOptions o=OptionsFor(si, "Si 2x2x2 shifted MP (k=±¼)");
    o.densityEcut=20.0; o.imposeSymmetry=true; o.kShift=rvec3_t(0.5,0.5,0.5);
    SCFParams par=TightGates(60);
    EnvOverrides(o, par);
    GpwReport report("Si "+o.label, par.Verbose);
    qchem::SolidCalculation calc(lat, MakeBasisSR(*si.cell), o, par);
    auto R=calc.Result();
    ASSERT_TRUE(R) << "the shifted mesh must converge, not merely stop: " << Why(R);
    EXPECT_NEAR(R->TotalCharge(), 8.0, 1e-6);
    // vs CP2K's OWN shifted 2x2x2 deck, re-measured through scripts/bench 2026-08-19 at -7.867436530436260.
    // Measured here: -7.868473428 (16 iterations, drho 1.0e-9) -- 1.04 mHa below, and the tolerance is the
    // historical 3 mHa.  Anything near -3.7 means the quarter-integer screen defect is back.
    EXPECT_NEAR(R->Energy(), -7.86744, 3e-3) << "GPW 2x2x2 shifted MP (CP2K default) vs -7.86744";
}


//  STEP 4's MIXED-MESH acceptance, at the smallest genuinely mixed mesh: N=(3,1,1) has ONE TRIM point
//  (Γ -- 2·1/3 is no reciprocal-lattice vector) beside two complex blocks (k=±1/3), so ONE composite
//  carries both child scalars through the full SCF -- the case Γ-only cannot reach.  (A 3×3×3 run is
//  the same code path 27 blocks wide; this keeps the gate's wall-time at ~3 Γ runs.)
TEST(GPW_Si, k311_eqCplx)
{
    const Material si=qchem::Materials::Get("Si_diamond");
    const Lattice_3D lat=LatticeOf(si, ivec3_t(3,1,1));
    ExpectRealComplexTwins(lat, si, "Si (3,1,1) mixed-mesh real-vs-complex");
}


//  MOM ON A MIXED MESH -- the END-TO-END gate for the R2.21 state split (flagged at that item's
//  landing, added here at the merge).  Until R2.21, a real TRIM block had nowhere to put its
//  mat_t<double> MOM reference (the references were typed by the RUN), so switching MOM on made a
//  flipped run THROW mid-SCF -- after the basis, seed and first Fock were already paid for.  Now each
//  block's reference lives in the shared OccupationState under its OWN scalar, so the real child fills
//  under a genuine OccupationPolicy<double> and the run keeps one ledger.
//
//  The gate drives the case that could not previously run AT ALL: (3,1,1) -- one real Γ block beside
//  two complex ±⅓ blocks -- with MOM armed early (MOMStartIter=2, so the reference is captured and
//  consulted well inside the run), real vs forceComplex.  Equal converged energy = the real block's
//  reference was captured, scored and applied exactly as its complex twin's.
TEST(GPW_Si, k311_Uni_MOM_eqCplx)
{
    const Material si=qchem::Materials::Get("Si_diamond");
    const Lattice_3D lat=LatticeOf(si, ivec3_t(3,1,1));

    SCFParams par;
    par.NMaxIter=60; par.MinΔρ=1e-3; par.MinΔE=1e-6;
    par.MinΔFD=1e30; par.MinVirial=1e30; par.MinFD=1e30; par.StartingRelaxRo=0.3; par.MergeTol=1e-4;
    par.UseMOM=true; par.MOMStartIter=2;              // armed EARLY: the reference must be live in-run

    qchem::SolidCalculation on (lat, MakeBasisSR(*si.cell),
                                {.Nelec=8, .species={{"Si",4}}, .densityEcut=20.0}, par);
    qchem::SolidCalculation off(lat, MakeBasisSR(*si.cell),
                                {.Nelec=8, .species={{"Si",4}}, .densityEcut=20.0, .forceComplex=true}, par);
    auto ron=on.Result(), roff=off.Result();
    ASSERT_TRUE(ron)  << "MOM on a mixed real/complex mesh must converge (R2.21): "
                      << (ron ? std::string() : ron.Error().details);
    ASSERT_TRUE(roff) << (roff ? std::string() : roff.Error().details);
    EXPECT_NEAR(ron->Energy(), roff->Energy(), 1e-9)
        << "a real TRIM block's MOM reference must behave exactly like its complex twin's";
    EXPECT_NEAR(ron->TotalCharge(), roff->TotalCharge(), 1e-10);
}
