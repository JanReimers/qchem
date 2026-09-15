// File GPW_SCF_UT.C  The GPW self-consistent total energy: the first periodic SCF on GAUSSIAN orbitals.
//
// GPW (increments 1-2) already satisfies every plane-wave Kohn-Sham concept EXCEPT the external potential:
//   - kinetic  -> Kinetic<dcmplx> calls bs->Kinetic()         (GPW: lattice-sum <p^2>)              [inc 1]
//   - Hartree  -> PW_Hartree casts bs to Orbital_DFT_IBS<dcmplx> + cd to FourierDensity (GPW: collocation tensors)[inc 2]
//   - XC       -> Vxc_Quadrature, same casts + the fit-basis grid                                          [inc 2]
//   - ion-ion  -> IonIon<dcmplx> (Ewald from Zion)                                                   [structure]
// The one gap was the external pseudopotential: the plane-wave PW_Pseudo needs G-space form factors, which
// Gaussians cannot supply.  This increment closes it: GPW_IBS realises Integrals_Pseudo<dcmplx> by REAL-SPACE
// mesh quadrature of the pseudopotential against its Gaussians (the SAME qcMesh machinery the molecular
// PP_Local/PP_NonLocal terms use).  So the ENTIRE Ham_PW_DFT drives a GPW basis verbatim -- Gaussian orbitals,
// plane-wave-style Hartree/XC by collocation -- through the real framework cSCFIterator.
//
// Validation (mirrors L_PP's finite==lattice cross-check, lifted to a full SCF):
//   (1) A single Si pseudo-atom in a LARGE cubic box, run through the GPW SCF, reproduces the FINITE molecular
//       density-fit DFT energy (the qchem::Calculation "sipp" + GTH-LDA pseudo-atom) to grid-cutoff tolerance.
//       Same basis, same PP, same functional; the only differences (density-fit vs collocation Hartree, Becke
//       vs uniform-grid XC, periodic vs open boundary) vanish as the box grows + the grid resolves -> the
//       electronic energies converge.  This is the tight correctness gate.
//   (2) A real material: crystalline silicon (diamond) at Gamma converges, conserves charge (8 valence e-),
//       and lands a reproducible total energy (a "did-E-move" regression anchor, per doc/GPWPlan.md section 5).
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
TEST(GPW_SCF, SiliconMultiKPlumbing)
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
TEST(GPW_SCF, SR_2x2x2GammaCentred_vs_CP2K)
{
    const Material si=qchem::Materials::Get("Si_diamond");
    const Lattice_3D lat=LatticeOf(si, ivec3_t(2,2,2));
    SolidCalcOptions o=OptionsFor(si, "Si 2x2x2 Gamma-centred");
    o.densityEcut=20.0; o.imposeSymmetry=true;
    SCFParams par=TightGates(60);
    EnvOverrides(o, par);
    GpwReport report("Si "+o.label, par.Verbose);
    qchem::SolidCalculation calc(lat, MakeBasisSR(*si.cell), o, par);
    auto R=calc.Result();
    ASSERT_TRUE(R) << Why(R);
    EXPECT_NEAR(R->TotalCharge(), 8.0, 1e-6);
    EXPECT_NEAR(R->Energy(), -7.77846, 3e-3) << "GPW 2x2x2 Gamma-centred vs CP2K same-mesh -7.77846";
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
TEST(GPW_SCF, SR_2x2x2ShiftedMP_vs_CP2K)
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
TEST(GPW_SCF, TermTranslationInvariance)   // RE-ENABLED 2026-09-15: 0.3 s, invariant to 1e-10 -- it had no live reason to be parked
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

// (1) THE GAMMA ANCHOR == THE CP2K ENERGY GATE.  SR basis, Rcut=2a (every term translation-invariant and
// screened-complete), densityEcut=20 (FFT N=32): reproduces the CP2K FCC-Si Gamma GPW reference (SIPP_SR /
// GTH-PADE-q4 / LDA_X+VWN5) Etot=-7.11506 to the N=32 grid gap (~0.4 mHa; densityEcut>=30 -> -7.11505 exact).
// NOTE (analytic path): the old fast Rcut=0 anchor (-8.2476) is GONE -- the analytic collocation is always
// screened-complete Bloch, so home-only 1E matrices would MIX SCHEMES (Tr(D S_home)=8 while the grid density
// integrates the Bloch trace -- the forbidden inconsistency; doc/GPWPlan.md durable pins).  SR keeps the Bloch
// overlap cleanly PD at 2a.  Energy-gated at the density-fit floor (minDE=1e-6, minDrho relaxed to 1e-3).
TEST(GPW_SCF, SiliconGammaConverges)
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

// The GPW run-report SCHEMA CHECK (RunReportPlan step 3).  Under an open run report the facade emits the
// `basis` (conditioning pre-flight) and `grids` (the ladder) sections itself during construction, and
// MakeIrrepWFs fills basis.perIrrep (per-Bloch-block conditioning) via the cursor.  Only the SETUP matters
// here, so a couple of iterations is plenty.
TEST(GPW_SCF, GridsReportSchema)
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
TEST(GPW_SCF, RealTRIMBlocksRunRealInReport)
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

// (2) THE TIGHT CROSS-CHECK: the isolated Si pseudo-atom in a box vs the finite molecular DFT on the SAME
// SIPP basis + GTH-LDA PP.  With the G-space local PP the GPW total is box-independent and reproduces the
// finite SIPP energy (-3.74 vs -3.759) to grid tolerance -- the doc/GPWPlan sec 3.4 correctness gate.
// NOTE: at Gamma the atom has NO point group, so its half-filled 3p shell is degenerate -- the ENERGY converges
// (-3.736, grid-stable) but the density rotates freely within that degenerate shell, so |Delta rho| never
// reaches the tolerance (not a bug: integer occupation of a degenerate open shell).  A dcmplx GDM/Ladder
// energy-minimiser would converge it (today GDM/Ladder are <double>-only); the crystal above sidesteps it with
// a gap.  So this pins the CONVERGED ENERGY + charge as a did-E-move anchor, without a Converged() guard.
TEST(GPW_SCF, SiPseudoAtomInBoxMatchesFinite)
{
    // Basis-MATCHED reference: the SAME SIPP Gaussian basis + GTH-LDA PP as a finite molecule (density-fit
    // Hartree, Becke XC).  This is the tight cross-check: GPW-in-box == finite molecular DFT (doc/GPWPlan sec 3.4).
    Molecule si; si.Insert(new Atom(14, 0.0, {0,0,0}));
    Calculation cSipp(si, {.basis = "sipp", .pseudopotential = true});
    const double Esipp=cSipp.Energy();
    // Physical oracle (near-complete Slater/High, for context -- a different basis, not the GPW-correctness gate).
    AtomCalculation cHi(14, 14-4, {.type=AtomType::Slater, .accuracy=BasisSetAccuracy::High, .pseudopotential=true});
    std::cout << "[Si finite] sipp="<<Esipp<<"  Slater/High="<<cHi.Energy()<<std::endl;

    // Box a=16 (was 11): the analytic collocation always includes the screened cross-cell pair products, so
    // the box must be large enough that they are negligible for the finite-molecule comparison (SIPP's most
    // diffuse alpha=0.06 pair prefactor: e^{-0.03 a^2} = 2.7e-2 at a=11 -- visible; 4.6e-4 at a=16).
    const Material box=qchem::Materials::Get("Si_box16");
    const Lattice_3D lat=LatticeOf(box);
    // XC route pinned UNIFORM (the gate's calibrated arrangement): under the Becke default this
    // DEGENERATE half-filled 3p atom exposed a real open question — the freely-rotating degenerate
    // density has orientation-DEPENDENT quadrature error on the fixed-axis Becke angular grid (V_xc is
    // nonlinear, so an anisotropic rho's error rotates with it), turning the energy-neutral rotation
    // into a ~Ha-scale E oscillation.  The SMEARED sibling (SmearingConvergesDegenerateShell) is fine
    // under Becke — fractional occupation restores the symmetric density.  Recorded in doc/GPWPlan1.md
    // (Becke remaining increments); this gate's PURPOSE is the PP box-independence check, so it keeps
    // its historical route.
    SolidCalcOptions o=OptionsFor(box, "Si atom-in-box");
    o.densityEcut=10.0; o.imposeSymmetry=true;
    o.seed=qchem::ChargeDensity::SeedStrategy::Uniform;
    o.images=BasisSet::Gaussian::CellImages::HomeCellOnly;    // the finite-molecule mode
    o.xcMesh.cellKind=qcMesh::UnitCellKind::Uniform;
    SCFParams par=TightGates(40);
    EnvOverrides(o, par);
    GpwReport report("Si "+o.label, par.Verbose);
    qchem::SolidCalculation calc(lat, MakeBasis(*box.cell), o, par);

    // Energy-converged; the density is degenerate at Gamma (see the note above), so these are the LAST
    // ITERATE's numbers by design -- the facade names them so, and this gate asks for exactly that.
    EXPECT_NEAR(calc.LastIterateCharge(), 4.0, 1e-6);         // 4 valence electrons (Zion=4), charge conserved
    // GPW-in-box (G-space local PP -> box-independent) reproduces the finite SIPP DFT energy to grid tolerance.
    EXPECT_NEAR(calc.LastIterateTerms().GetTotalEnergy(), Esipp, 5e-2) << "GPW-in-box total vs finite SIPP molecular DFT";
}

// (tier 4b, invariant) THE ζ=0 COLLAPSE: the TWO-CHANNEL machinery on a CLOSED SHELL must reproduce the
// unpolarized anchor.  Same gapped Si/Γamma cell + recipe as SmearingInertOnGap, but multiplicity=1 drives
// the polarized pipeline (the dcmplx composite WF under SpinGroup::Polarized, Crystal_EC(4,4), the spin-native XC term) with
// nUp=nDn=4 -- v^σ(ρ/2,ρ/2)=v^P(ρ) pointwise, so the total must land on the SAME −7.11506 anchor.  The
// periodic sibling of the molecular WaterPolarizedLDA-vs-LDA check; catches any polarized-path divergence
// (channel bookkeeping, shared-engine caching, the collocation memo screen) on known ground.
TEST(GPW_SCF, PolarizedSingletMatchesUnpolarizedSiGamma)
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
TEST(GPW_SCF, SharedFermiLevelLetsTheMomentRelax)
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
    // The moment relaxed away: this lands on the singlet anchor of PolarizedSingletMatchesUnpolarizedSiGamma.
    EXPECT_NEAR(shared.E, -7.11506, 3e-3)
        << "a shared mu must let a triplet-seeded closed-shell system fall back to the singlet";

    // CONTROL: the same run with separate per-channel counts cannot relax -- nUp-nDn=2 is conserved.
    const Arm held=TripletSi(false);
    EXPECT_NEAR(held.charge, 8.0, 1e-6);
    EXPECT_GT(held.E, shared.E + 1e-3)
        << "with two reservoirs the seeded multiplicity is a CONSTRAINT and must sit above the free minimum";
}

TEST(GPW_SCF, PolarizedSeedSingletMatchesUnpolarizedSiGamma)
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

// (tier 4b, gate b) O2 in a box, TRIPLET: the multi-electron polarized solid pipeline vs the finite
// molecular facade on the SAME sipp O-q6 basis + GTH PP (the spin sibling of SiPseudoAtomInBoxMatchesFinite,
// cross-anchored to the facade's spin-native triplet machinery -- doc/SymmetryUpgradePlan.md §4 tier 4b).
TEST(GPW_SCF, O2TripletInBoxMatchesFinite)
{
    const double d=2.282;   // O2 bond (au)
    Molecule o2mol;
    o2mol.Insert(new Atom(8, 0.0, {-d/2,0,0}));
    o2mol.Insert(new Atom(8, 0.0, { d/2,0,0}));
    Calculation cRef(o2mol, {.basis="sipp", .multiplicity=3, .pseudopotential=true});
    const double Eref=cRef.Energy();
    {
        auto E=cRef.EnergyTerms();
        std::cout << "[O2 finite] sipp GTH-q6 LSDA triplet="<<Eref
                  << "  (Ekin="<<E["Kinetic"]<<" Een="<<E["Een"]<<" Eee="<<E["Eee"]<<" Exc="<<E["Exc"]
                  << " Enn="<<E["Enn"]<<")"<<std::endl;
    }

    const Material box=qchem::Materials::Get("O2_box16");   // the same d=2.282 dimer, centred in a 16-bohr box
    const Lattice_3D lat=LatticeOf(box);
    SolidCalcOptions o=OptionsFor(box, "O2 in-box triplet");
    o.imposeSymmetry=true;   // V1.30: was the DEFAULT; now stated, because an imposition you did not ask for is invisible in the result
    o.multiplicity=3;                                  // S=1: nUp=7, nDown=5; densityEcut stays AUTO: O q6 is hard (alpha_max rules)
    o.images=BasisSet::Gaussian::CellImages::HomeCellOnly;
    GpwReport report("O "+o.label, false);
    qchem::SolidCalculation calc(lat, MakeBasis(*box.cell), o, TightGates(60));
    // No convergence guard (as before): the box's energy is what is compared, and it is compared as the
    // last iterate, which is what this gate has always read.
    EXPECT_NEAR(calc.LastIterateCharge(), 12.0, 1e-6);
    EXPECT_NEAR(calc.LastIterateTerms().GetTotalEnergy(), Eref, 5e-2) << "GPW-in-box triplet vs finite molecular LSDA triplet";
}

// (tier 4b, gate a) THE POLARIZED SOLID PIPELINE: Na pseudo-atom in a box, DOUBLET (doc/SymmetryUpgradePlan.md
// §4).  The minimal end-to-end TWO-CHANNEL GPW run: Na q1 GTH PP, S=1/2, moment 1 -- spin-resolved D through
// Crystal_EC(nUp=1,nDown=0), the dcmplx composite WF under SpinGroup::Polarized (two Bloch channels), and the spin-native Becke XC
// term.  Cross-anchored against the finite molecular facade doublet on the SAME
// valence basis + PP (the spin sibling of SiPseudoAtomInBoxMatchesFinite).
//
// SEED PIN (the 2026-08-04 root-cause campaign): this gate MUST seed from IonicSAD.  From the Uniform seed
// the lone ↑ electron converges to a GENUINE self-consistent excited basin 72 mHa above the minimum
// (diffuse 3s; DIIS honors it, GDM -- a local descender -- stays in it; total is box/grid/route-independent,
// so it looks "converged" by every health metric).  The functional itself was proven correct everywhere:
// fixed-density term probes (DISABLED_NaFixedDensityTermProbe) match analytic kinetic, the exact discrete
// G!=0 lattice-sum Hartree, and the ζ=1 Dirac/VWN values exactly, and E_GPW[D*]=−0.1420 at the independent
// radial same-basis oracle's minimizer (oracle E=−0.1416; complete-basis −0.1922).  A ONE-ELECTRON system
// is uniquely basin-fragile: no other electrons pull the density into the core basin (Na2, O2-triplet and
// F-doublet all escape Uniform fine).  IonicSAD lands in-basin: −0.1419, 4.8 mHa from the facade.
// SPIN-SAD (§10 increment B): the polarized run now seeds the TWO-CHANNEL PolarizedSeedCD, and Na's library
// pair is exact (1 e-: up=total, dn=0), so iteration 0 starts FULLY polarized -- same basin, same pin
// (-0.141933 to the digit), 21 iters vs the rho/2-collapse seed's 14 (the staggered start makes DIIS
// reorganize more; basin selection, not speed, is what the seed pin protects).
TEST(GPW_SCF, NaPseudoAtomInBoxDoublet)
{
    // Basis-MATCHED finite reference: the SAME valence_lowq_sr Na basis + GTH-LDA q1 PP through the
    // molecular facade, spin-native LSDA doublet (Ham_PP polarized: FittedVxcPol + FittedVcorrPol).
    // ppValence=1 overrides Na's GTH DEFAULT-valence entry, which is the SEMICORE q9 -- the q1 entry is
    // the valence_lowq bases' convention.  (The Slater/High atom path is NOT a good oracle here: Slater
    // functions fit the smooth nodeless pseudo-orbitals poorly -- same reason the Si gate anchors on the
    // basis-matched sipp facade run.)
    Molecule na; na.Insert(new Atom(11, 0.0, {0,0,0}));
    Calculation cRef(na, {.basis="valence_lowq_sr", .multiplicity=2, .pseudopotential=true, .ppValence=1});
    const double Eref=cRef.Energy();
    std::cout << "[Na finite] valence_lowq_sr LSDA doublet (q1)="<<Eref<<std::endl;

    // Box a=16 (the Si gate's size: cross-cell products of the most diffuse pair negligible).
    const Material box=qchem::Materials::Get("Na_box16");
    const Lattice_3D lat=LatticeOf(box);
    SolidCalcOptions o=OptionsFor(box, "Na atom-in-box doublet");
    o.imposeSymmetry=true;   // V1.30: was the DEFAULT; now stated, because an imposition you did not ask for is invisible in the result
    o.multiplicity=2;                                          // S=1/2: nUp=1, nDown=0
    o.images=BasisSet::Gaussian::CellImages::HomeCellOnly;   // the finite-molecule mode
    o.seed=qchem::ChargeDensity::SeedStrategy::IonicSAD;       // SEED PIN: Uniform has a stable wrong basin (header)
    GpwReport report("Na "+o.label, false);
    qchem::SolidCalculation calc(lat, MakeBasisLowQ(*box.cell, BasisSetData::VALENCE_LOWQ_SR), o, TightGates(40));
    auto R=calc.Result();
    ASSERT_TRUE(R) << "3s^1 is non-degenerate: Δρ converges -- " << Why(R);
    EXPECT_NEAR(R->TotalCharge(), 1.0, 1e-6);                // 1 valence electron (Zion=1), charge conserved
    // GPW-in-box (two-channel) reproduces the finite molecular LSDA doublet (measured 4.8 mHa; the gap is
    // the two stacks' fit/quadrature tech -- the facade Dunlap-fits J and fits v_xc, the GPW is fit-free).
    EXPECT_NEAR(R->Energy(), Eref, 2e-2) << "GPW-in-box doublet vs finite molecular LSDA doublet";
    EXPECT_NEAR(R->Energy(), -0.141933, 1e-4);               // did-E-move anchor (== the same-basis oracle -0.1416)
}

// (4b-i) FERMI SMEARING IS INERT ON A GAP (doc/GPWPlan1.md 4b, gate i).  The same gapped Si/Gamma anchor as
// SiliconGammaConverges, but with smearing kT=1e-3 Ha turned ON.  Si is a wide-gap insulator in this basis
// (ε_LUMO−ε_HOMO ≫ kT), so every f_i is essentially 0 or 1: the fractional occupations collapse to the
// aufbau integers, the Mermin −TS is negligible, and the total reproduces the CP2K reference −7.11506.  This
// is the regression that smearing must not perturb a system that does not need it (the T→0 / kT≪gap limit).
TEST(GPW_SCF, SmearingInertOnGap)
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

// (4b-iii) FERMI SMEARING CONVERGES A DEGENERATE OPEN SHELL (doc/GPWPlan1.md 4b, gate iii + the cure).  The
// Si pseudo-atom in a box has 4 valence electrons in a 3s²3p² configuration; at Gamma the atom has NO point
// group, so its three 3p orbitals are EXACTLY degenerate and half-filled.  Integer aufbau must pick 2 of the
// 6 p-states arbitrarily -> the density rotates freely within the degenerate shell and |Δρ| never converges
// (the documented behaviour of SiPseudoAtomInBoxMatchesFinite, iters=40/Δρ=0.08/"DENSITY-DEGENERATE", which
// pins only the energy).  Fermi smearing is the cure: μ lands in the degenerate manifold, each 3p orbital
// takes the SAME fractional occupation, the density is symmetric and STATIONARY, and the SCF converges Δρ
// (iters=24, Δρ=9e-7, "CONVERGED").  The Mermin −TS<0, so the total GetTotalEnergy() reported IS the free
// energy A=E−TS, which sits below the internal energy E -- the finite-T thermodynamic ordering (gate iii).
// kT MUST exceed the near-degenerate splitting to stabilise: kT=1e-2 converges, kT=1e-3 still slosh-rotates
// (measured) -- the honest recipe is "smear wider than the frontier splitting you are curing".
TEST(GPW_SCF, SmearingConvergesDegenerateShell)
{
    Molecule si; si.Insert(new Atom(14, 0.0, {0,0,0}));
    Calculation cSipp(si, {.basis = "sipp", .pseudopotential = true});
    const double Esipp=cSipp.Energy();   // finite SIPP molecular DFT (symmetric occupation): -3.759

    const Material box=qchem::Materials::Get("Si_box16");
    const Lattice_3D lat=LatticeOf(box);
    SolidCalcOptions o=OptionsFor(box, "Si atom-in-box +smear");
    o.densityEcut=10.0; o.imposeSymmetry=true;
    o.seed=qchem::ChargeDensity::SeedStrategy::Uniform;
    o.images=BasisSet::Gaussian::CellImages::HomeCellOnly;
    SCFParams par=TightGates(60);
    par.SmearingkT=1e-2;
    EnvOverrides(o, par);
    GpwReport report("Si "+o.label, par.Verbose);
    qchem::SolidCalculation calc(lat, MakeBasis(*box.cell), o, par);
    auto R=calc.Result();
    ASSERT_TRUE(R) << "Fermi smearing should converge Δρ where integer aufbau cannot (degenerate 3p): " << Why(R);
    EXPECT_NEAR(R->TotalCharge(), 4.0, 1e-6);
    const qchem::EnergyBreakdown E=R->EnergyTerms();
    // did-E-move anchor: the converged free energy A=E−TS at kT=1e-2 (internal E≈-3.744; the ~38 mHa gap to A
    // is the 3p-shell entropy −TS at this kT, which lowers A below E and below Esipp).  (Re-pinned when the
    // field-sharpness density rule landed -- doc/GPWPlan1.md 4b: the sharper XC grid moved it -3.779 -> -3.783.)
    EXPECT_NEAR(E.GetTotalEnergy(), -3.78260, 3e-3);
    EXPECT_LT(E["MinusTS"], 0.0);                               // −TS<0 => A=GetTotalEnergy() sits below internal E (gate iii)
    EXPECT_NEAR(E.GetTotalEnergy()-E["MinusTS"], Esipp, 3e-2) << "internal E=A−(−TS) vs finite SIPP molecular DFT";
}

// (item 2a) THE MOTIVATION: integer aufbau CANNOT converge the degenerate 3p.  With smearing OFF, aufbau must
// place the lone 3p electron in ONE of the three degenerate p orbitals -- an arbitrary, symmetry-broken pick.
// The TOTAL ENERGY settles (|ΔE/E|~1e-13, the gap column ~0 => no frontier gap, the metallic signature) but
// the DENSITY rotates freely within the degenerate manifold, so |Δρ| floors well above tolerance and never
// converges.  This is the honest reason the smearing/annealing path below exists (mirrors the documented
// SiPseudoAtomInBoxMatchesFinite degenerate-shell behaviour, now for a periodic lattice).
TEST(GPW_SCF, AlFCCDegenerateShellAufbauStalls)
{
    const Material al=qchem::Materials::Get("Al_fcc");        // Al (Zion=3): 3s^2 3p^1
    const Lattice_3D lat=LatticeOf(al);
    SolidCalcOptions o=AlOptions(al);
    SCFParams par=AlGates(40); par.SmearingkT=0.0; par.Verbose=true;   // aufbau (no smearing)
    GpwReport report("Al "+o.label, par.Verbose);
    qchem::SolidCalculation calc(lat, MakeBasisLowQ(*al.cell, BasisSetData::VALENCE_LOWQ_SR), o, par);

    EXPECT_FALSE(calc.Result()) << "integer aufbau cannot converge Δρ of a partially-filled degenerate 3p shell";
    EXPECT_NEAR(calc.LastIterateCharge(), 3.0, 1e-6);       // charge is still conserved (3 valence e-)
    EXPECT_NEAR(calc.LastIterateTerms()["MinusTS"], 0.0, 1e-12);   // no smearing => no entropy term
}

// (item 2b) THE CURE + ANNEALING (doc/GPWPlan1.md item 2): a DESCENDING kT schedule (0.02 -> 0.01 -> 0.005 Ha),
// re-seeding each stage from the previous converged density.  Fermi smearing puts μ in the degenerate manifold
// so each 3p orbital takes the SAME fractional occupation -- the density is cubic-symmetric and STATIONARY, and
// every stage converges Δρ where aufbau cannot.  Annealing makes the cold end cheap: the kT=0.005 cold-START
// takes ~44 iters, but re-seeded from kT=0.01 it converges in ~17.  The Mermin −TS<0, so GetTotalEnergy() is
// the free energy A=E−TS (below the internal E); the INTERNAL energy E=A−(−TS) is kT-INDEPENDENT to ~1e-7
// across all three stages (−1.92115) -- the physical T→0 answer, and a strong self-consistency check on the
// smearing thermodynamics (gate iii).
TEST(GPW_SCF, AlFCCAnnealedMetal)
{
    const Material al=qchem::Materials::Get("Al_fcc");        // Al (Zion=3): 3s^2 3p^1
    const Lattice_3D lat=LatticeOf(al);
    SolidCalcOptions o=AlOptions(al);
    // Accelerator = plain DIIS (NOT the DIIS->GDM Ladder).  MEASURED 2026-07-28: the Ladder's GDM tail rung is
    // INCOMPATIBLE with Fermi smearing here.  GDM builds its geodesic DIRECTION from the fixed-occupation
    // electronic gradient [F,D], but line-searches the FREE energy A=E−TS with the occupations Fermi-refilled
    // per trial (SCFIterator DirectMinStep).  Under fractional occupation those disagree: A is stationary where
    // the smeared gradient (not [F,D]) is zero, so at the DIIS fixed point [F,D]~4e-3≠0, and GDM's first step is
    // a persistent FALLBACK (GPW_GDMTRACE: best(Et−Ecur)=+0.42, all 12 backtracks uphill) -> occupations flip
    // (cfg `*`), the run never recovers.  This is NOT a grid / non-variationality fault (DIIS converges the SAME
    // E[ρ] cleanly to −1.97521); GDM just needs the occupation-response term to tail-polish a smeared free
    // energy.  So DIIS here; GDM+smearing is a captured follow-up (doc/GPWPlan1.md item 2).
    std::vector<qchem::SCFStage> schedule;
    for (double kT : {0.02, 0.01, 0.005})
    {
        SCFParams par=AlGates(); par.SmearingkT=kT;
        schedule.push_back({par, qchem::SCFAccelerators::Type::DIIS});
    }
    GpwReport report("Al "+o.label, false);
    qchem::SolidCalculation calc(lat, MakeBasisLowQ(*al.cell, BasisSetData::VALENCE_LOWQ_SR), o, schedule);
    auto R=calc.Result();
    ASSERT_TRUE(R) << "Fermi-smearing annealing converges Δρ where integer aufbau cannot (degenerate 3p): " << Why(R);
    EXPECT_NEAR(R->TotalCharge(), 3.0, 1e-6);
    const qchem::EnergyBreakdown E=R->EnergyTerms();
    EXPECT_LT(E["MinusTS"], 0.0);                              // −TS<0 => A=GetTotalEnergy() sits below internal E (gate iii)
    // did-E-move anchors at the coldest stage (kT=0.005): the free energy A and the kT-independent internal E.
    EXPECT_NEAR(E.GetTotalEnergy(), -1.934665, 2e-3);                          // A = E − TS at kT=0.005
    EXPECT_NEAR(E.GetTotalEnergy()-E["MinusTS"], -1.921148, 2e-3);              // internal E (T→0 physical value)
}

// (item 3) GLOBAL μ ACROSS k-BLOCKS -- the true metal fill.  FCC Al on a 2×2×2 Γ-centred Bloch mesh (8
// k-blocks) with ONE chemical potential across the whole BZ (Crystal_EC global mode + the composite cross-k
// Fermi fill): charge SLOSHES between k-points under a single μ instead of each k being pinned to a fixed
// per-block count.  This is the structural step a metal needs -- the partially-filled 3p band disperses with
// k, so no per-k integer (or per-k Fermi) occupation is right; the physical occupation is set by where each
// k's bands sit relative to the ONE Fermi level.  WHAT THIS GATES: (a) charge is conserved as the BZ-weighted
// Σ_k w_k n_k = 3 (the weight-consistency guard -- the μ constraint uses the SAME w_k the density applies);
// (b) the single μ CONVERGES the mesh where per-block filling cannot (measured: AL_GLOBAL=0 at 2×2×2 forces 3
// e⁻ at every k and lands non-converged garbage A≈-0.46, vs the global μ's converged -2.117 -- the charge
// MUST redistribute between k-points); (c) k-sampling lowers the energy vs Γ-only (-1.92 → -2.12, real
// dispersion).  Reduces EXACTLY to the per-block Fermi at a single k (verified: global≡per-block to 1e-12 at Γ).
TEST(GPW_SCF, AlFCCMetalGlobalMu)
{
    const Material al=qchem::Materials::Get("Al_fcc");        // Al (Zion=3): 3s^2 3p^1
    const Lattice_3D lat=LatticeOf(al, ivec3_t(2,2,2));       // 8-point Γ-centred Bloch mesh (weights sum to 1)
    SolidCalcOptions o=AlOptions(al, "Al FCC 2x2x2 global mu");
    o.imposeSymmetry=false;                    // the FULL mesh: the IBZ-folded sibling below must reproduce it
    o.globalFermi=true;                        // ONE μ across the BZ (the metal)
    SCFParams par=AlGates(); par.SmearingkT=0.01;
    GpwReport report("Al "+o.label, false);
    qchem::SolidCalculation calc(lat, MakeBasisLowQ(*al.cell, BasisSetData::VALENCE_LOWQ_SR), o, par);
    auto R=calc.Result();
    ASSERT_TRUE(R) << "one μ across the BZ converges the dispersive metal (per-block filling cannot): " << Why(R);
    EXPECT_NEAR(R->TotalCharge(), 3.0, 1e-6);  // BZ-weighted Σ_k w_k n_k = 3 (weight-consistency guard)
    const qchem::EnergyBreakdown E=R->EnergyTerms();
    EXPECT_LT(E["MinusTS"], 0.0);              // −TS<0 => A=GetTotalEnergy() is the free energy (gate iii)
    std::cout<<"[Al global-μ full-mesh] A="<<E.GetTotalEnergy()<<std::endl;   // the IBZ pin's route-matched partner
    EXPECT_NEAR(E.GetTotalEnergy(), -2.11681, 3e-3);   // did-E-move anchor (2×2×2 global-μ free energy A)
    EXPECT_LT(E.GetTotalEnergy(), -1.95);      // dispersion: well below the Γ-only -1.92 (k-sampling binds)
}

// (item 3, IBZ/k-star) FOLDING IS EXACT.  The 8-point 2×2×2 Γ-mesh folds to 3 irreducible k-points under the
// cubic point group, and the density is STAR-AVERAGED consistently: the G-space Hartree via the reciprocal ops
// the basis exposes (GetReciprocalPointOps, ctor-injected into the composite density), and the real-space XC
// raster via the DIRECT ops (cFIT_SF_ABS::SymmetrizeRaster, ctor-injected into the Vxc fit basis -- so XC stays
// on the non-negative ρ_DM raster).  So the reduced run reproduces the full-mesh AlFCCMetalGlobalMu free energy
// to grid/SCF tolerance (measured ~6e-8) with fewer k-points -- the IBZ payoff, done exactly (doc/GPWPlan1 item 3).
TEST(GPW_SCF, AlFCCMetalIBZExact)
{
    const Material al=qchem::Materials::Get("Al_fcc");
    const Lattice_3D lat=LatticeOf(al, ivec3_t(2,2,2));
    SolidCalcOptions o=AlOptions(al, "Al FCC 2x2x2 IBZ");
    o.globalFermi=true; o.imposeSymmetry=true;       // fold to the irreducible wedge AND star-average the density
    SCFParams par=AlGates(); par.SmearingkT=0.01;
    GpwReport report("Al "+o.label, false);
    qchem::SolidCalculation calc(lat, MakeBasisLowQ(*al.cell, BasisSetData::VALENCE_LOWQ_SR), o, par);
    auto R=calc.Result();
    ASSERT_TRUE(R) << Why(R);
    EXPECT_NEAR(R->TotalCharge(), 3.0, 1e-6);
    // Re-anchored 2026-08-01: the 0i custom V_loc G-ball (harmonic routing + custom top level) moved the
    // long-PP block by 1.6e-4 on Al's coarse grids -- full mesh AND reduced shift TOGETHER (folding stays
    // exact).  Old kappa-sweep anchor: -2.116812.
    // Re-pinned 2026-08-08 ON THE UNIFORM ROUTE (V2.4, arming the V1.26 cost selector): Al is soft
    // (alpha_max=4), so Auto now costs the two grids and picks uniform -- 4,096 points against Becke's
    // 18,000.  The value it lands on, -2.1169707, is EXACTLY the uniform-route number the previous
    // comment here already recorded, so this re-pin introduces no new quantity; it switches which of two
    // long-known values is the default.
    //   Which is right?  Measured, not assumed (GPW_SCF.DISABLED_GridRouteAB_AlFCC): against a fine Becke
    //   reference the uniform route is CLOSER on both scores -- ||drho||_1 6.18e-4 vs the production
    //   Becke mesh's 8.95e-4, and dEtot -1.34e-5 vs +1.93e-4.  The production Becke mesh's own residual
    //   is the angular error V2.6 measured on this system (Al is the worst case there).  And the uniform
    //   route is converged: refining its cutoff 4x moves ||drho||_1 by 1 part in 6000.
    // Previous anchor, on the Becke route: -2.1174805.
    EXPECT_NEAR(R->Energy(), -2.1169707, 1e-4)           // == the full 8-k-point mesh: IBZ symmetrization is EXACT
        << "IBZ-reduced must reproduce the full-mesh free energy (AlFCCMetalGlobalMu prints the same value)";
}

// (item 5, IBZ) NON-SYMMORPHIC -- diamond Si (FCC lattice + a 2-atom basis at (0,0,0),(¼,¼,¼)) is space group
// Fd-3m: NON-symmorphic (the two sublattices are related by a glide, τ=(¼,¼,¼)≠0).  The k-FOLD reaches the full
// Oh (3 irreducible k-points, same as FCC Al on this lattice): the τ=0 Td subgroup + time reversal (k→−k) already
// supplies the inversion Td lacks.  The DENSITY star-average now carries the glide τ: the G-space Hartree via the
// e^{+2πi(Um)·τ} phase (SymmetrizeGMap over SpaceGroup::ReciprocalOps) and the real-space XC raster via the exact
// FFT fractional shift ρ(W·x+τ) (SymmetrizeRaster over SpaceGroup::DirectOps).  So the IBZ-reduced total now
// reproduces the full-mesh Γ-centred 2×2×2 exactly (was −8.259, ~0.48 Ha off, under the old τ=0 W-only guard).
TEST(GPW_SCF, SiDiamondIBZ_NonSymmorphic)
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

// (T3.2, doc/SymmetryUpgradePlan.md §6b) STREAM FOLD through-SCF A/B on an IMPOSED Γ-only run: the factory
// arms route (b) on the shared molecular evaluator (reduced stream build + replay, rep-transform h), and the
// existing SymmetrizeGMap/SymmetrizeRaster sites complete the group-average.  GPW_STREAM_FOLD=0/1 toggles the
// fold between two otherwise identical runs (read fresh in the factory, so one process can A/B).  The two
// totals must agree to the band-limit class (§8 through-SCF tier -- the production 5-smooth grid is NOT
// τ-commensurate, so reduced+P and full+P are two equally valid quadratures of the same density), and the
// folded run's [stream cache] line must show the reduced build (repPairs << pairs).
TEST(GPW_SCF, StreamFoldImposedGamma_SiDiamond)
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

// (T3.5, doc/SymmetryUpgradePlan.md §6b) THE ARMING GATE: the case that FORCED the 2026-08-03 default-on
// retraction, re-run as a standing A/B.  The Si pseudo-atom in a box at Γ is a DEGENERATE OPEN SHELL --
// 3s²3p² with three exactly-degenerate p orbitals holding two electrons, at integer aufbau (no smearing).
// Its density therefore rotates freely inside the degenerate manifold and NEVER converges Δρ (the documented
// SiPseudoAtomInBoxMatchesFinite behaviour; this gate pins ENERGY, like that one, and does not assert
// convergence).  Under the OLD fold the reduced replay SAMPLED each pair orbit's representative D element,
// which asserts a symmetric D that this run breaks permanently: the armed run flipped out of the benign
// rotating-ρ mode into charge-transfer sloshing, ~0.26 Ha off, and default-on was withdrawn.  With the
// replay reading the ORBIT-PROJECTED D (T3.5) the folded and unfolded imposed runs solve the same equations
// -- P ρ_red[D] = ρ[P D] = P ρ_full[D] for ANY iterate -- so the two arms must now agree to the band-limit
// class ON EXACTLY THIS CELL.  That agreement is the whole licence for arming the fold by default; if this
// gate ever reopens the 0.26 Ha gap, the default goes back to opt-in.
TEST(GPW_SCF, StreamFoldOpenShellMatchesUnfolded_SiAtomInBox)
{
    const Material box=qchem::Materials::Get("Si_box16");     // Pm-3m box, 48 ops; the atom sits on the cube centre
    const Lattice_3D lat=LatticeOf(box);                      // Γ-only: the T3.2 arming condition
    SolidCalcOptions o=OptionsFor(box, "Si atom-in-box Γ open-shell fold A/B");
    o.densityEcut=10.0; o.imposeSymmetry=true;
    o.images=BasisSet::Gaussian::CellImages::HomeCellOnly;    // the finite-molecule mode of the parent gate
    o.xcMesh.cellKind=qcMesh::UnitCellKind::Uniform;            // ditto: the rotating degenerate density and
    o.seed=qchem::ChargeDensity::SeedStrategy::Uniform;         //   a fixed-axis Becke grid do not mix
    SCFParams par=TightGates(40); par.SmearingkT=0.0;           // INTEGER AUFBAU -- the symmetry-broken D
    struct Arm { double charge, E; };
    auto arm=[&](const char* fold) -> Arm                       // last iterates: this gate pins ENERGY, not convergence
    {
        setenv("GPW_STREAM_FOLD",fold,1);  qchem::ReresolveRunPolicy();
        GpwReport report("Si "+o.label, false);
        qchem::SolidCalculation calc(lat, MakeBasis(*box.cell), o, par);
        return {calc.LastIterateCharge(), calc.LastIterateTerms().GetTotalEnergy()};
    };
    const Arm R0=arm("0"), R1=arm("1");
    unsetenv("GPW_STREAM_FOLD");  qchem::ReresolveRunPolicy();

    std::cout << "[open-shell fold A/B] E(full)=" << R0.E << "  E(folded)=" << R1.E << "  dE=" << R1.E-R0.E << std::endl;
    EXPECT_NEAR(R1.charge, 4.0, 1e-6);
    EXPECT_NEAR(R0.charge, 4.0, 1e-6);
    EXPECT_NEAR(R1.E, R0.E, 1e-3)
        << "a DEGENERATE OPEN SHELL must not care whether the streams are folded (the retracted 0.26 Ha)";
}

// ===== EXPERIMENTAL (scratch): global-μ across k-blocks (item 3 inc 3) =====
// AL_KGRID=n (mesh nxnxn), AL_GLOBAL=0/1 (per-block vs global μ), AL_KT, AL_NMAX.
// I2 (plan §6a fit/grid SEPARATION): the (Delta, uniform) cross cell.  The SAME material through the
// PLANE-WAVE fit (band-limited v_xc on the FFT raster) and through the DELTA fit on the uniform cell
// mesh -- the two v_xc representations must agree to the route-difference class (band-limiting +
// raster-vs-midpoint-mesh quadrature), the same class the Becke-vs-uniform gate measures (~1e-4 Exc).
TEST(GPW_SCF, DeltaFitUniformGridMatchesPWFit_SiGamma)
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
// unpolarized PW-fit answer, exactly as PolarizedSingletMatchesUnpolarizedSiGamma pins it on the Becke route.
TEST(GPW_SCF, PolarizedSingletMatchesUnpolarized_PWFitRaster)
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

// THE W1 GATE (doc/SymmetryUpgradePlan.md §6a): Becke XC under IBZ.  BeckeFit_IBS group-averages its
// mesh INVARIANT and star-averages rho every iteration (the SymmetrizeRaster hook, exact orbit-mean
// projector); on the converged SYMMETRIC density the invariant mesh integrates identically to the
// single-orientation mesh (Q_inv(f)==Q(f) for a symmetric f), so the IBZ run must reproduce the
// full-mesh Becke run to the IBZ class -- the Becke sibling of SiDiamondIBZ_NonSymmorphic, with the
// non-symmorphic glide tau exercised through the torus fold + MakeInvariant.  COARSE explicit Becke
// recipe on BOTH arms: the gate compares like against like, so grid quality cancels (and the
// imposed arm's group-average mesh growth stays affordable).
TEST(GPW_SCF, BeckeXC_IBZ_SiDiamond)
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

// (item 4) THE HONEST METAL: FCC Na, a real half-filled-band Fermi surface.  Zion=1 (3s^1) => ONE valence
// electron per cell, so the single conduction band is HALF-FILLED and μ cuts THROUGH it -- a genuine Fermi
// surface (unlike Al's degenerate-3p at Γ).  SHIFTED Monkhorst-Pack 2×2×2 (kShift=½ => k at ±¼, CP2K's default
// -- avoids the high-symmetry Γ, samples the Fermi surface evenly) + global μ + Fermi smearing.  MEASURED: μ
// lands mid-band, the 2 k-points inside the Fermi surface fill (n_k=2.0) while the 6 on it smear FRACTIONALLY
// (n_k=0.67, ε≈μ, f=1/(1+e^{0.7})=0.33/spin) -- the textbook smeared Fermi surface, charge Σ_k w_k n_k = 1
// exactly, converged in ~26 iters.  BASIS: VALENCE_LOWQ_SR2 Na at the REAL FCC-Na density (a=10 au, matched to
// Na's atomic volume) -- SR2 drops the diffuse Na s 0.0857 + p 0.05, so the Bloch overlap is well-conditioned
// at the correct lattice constant (cond~38; SR needs an unphysical a=12 to condition).  Remaining tension: SR2
// is a MINIMAL 6-function basis (no diffuse 3s), so a fuller metallic Na basis (valgen, step 1) + IBZ/mesh-
// convergence is the accuracy follow-up; this gate validates the machinery + Fermi surface, not a cohesive E.
TEST(GPW_SCF, NaFCCMetalGlobalMu)
{
    const Material na=qchem::Materials::Get("Na_fcc");        // FCC Na at Na's atomic density; 3s^1 => half-filled band
    const Lattice_3D lat=LatticeOf(na, ivec3_t(2,2,2));
    SolidCalcOptions o=OptionsFor(na, "Na FCC metal");
    o.imposeSymmetry=true;   // V1.30: was the DEFAULT; now stated, because an imposition you did not ask for is invisible in the result
    o.globalFermi=true;                        // ONE μ across the BZ
    o.kShift=rvec3_t(0.5,0.5,0.5);             // shifted Monkhorst-Pack (k at ±¼)
    o.seed=qchem::ChargeDensity::SeedStrategy::Uniform;
    SCFParams par=Gates(60, 1e-5, 1e30); par.SmearingkT=0.01;
    GpwReport report("Na "+o.label, false);
    qchem::SolidCalculation calc(lat, MakeBasisLowQ(*na.cell, BasisSetData::VALENCE_LOWQ_SR2), o, par);
    auto R=calc.Result();
    ASSERT_TRUE(R) << "global μ + smearing converges the half-filled-band metal: " << Why(R);
    EXPECT_NEAR(R->TotalCharge(), 1.0, 1e-6);  // one valence electron, BZ-weighted Σ_k w_k n_k = 1
    EXPECT_LT(R->EnergyTerms()["MinusTS"], -1e-4);   // −TS<0 AND non-trivial: the Fermi surface IS fractionally filled
    EXPECT_NEAR(R->Energy(), 0.045543, 3e-3);  // did-E-move anchor (free energy A at kT=0.01)
}

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
TEST(GPW_SCF, NaFRocksaltGamma)   // RE-ENABLED 2026-09-15: 15 s, converged in 23 iterations -- the GPW x NaF row's first standing anchor
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
TEST(GPW_SCF, DISABLED_NaFGridContinuation)
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

//================================================================================================
//  THE BECKE XC GATE (doc/GPWPlan1.md "Becke XC grid").  The atom-centred periodic Becke XC
//  quadrature must reproduce the uniform-multigrid XC on a CONDITIONED basis -- same E_xc and same
//  V_xc matrix to grid tolerance -- before it can become the default for diffuse bases.  Converge
//  Si/Gamma on the standard uniform route (the SiliconGammaConverges recipe), then evaluate BOTH
//  XC term pairs (Dirac + VWN5) on the SAME converged density:
//    uniform -- the PAIR quadrature on the Vxc fit basis's FFT grid (the raw-collocation route);
//    Becke   -- Vxc_Quadrature: rho(r) analytic at the atom-centred points, MatrixOverlap matrix.
//  Angular rule: GaussLegendre (machine-exact algebraic degree at any L -- the audited Lebedev
//  tables stop at L=11; see the Mesh_AngularDegree tests).
//================================================================================================
TEST(GPW_SCF, BeckeXCMatchesUniformXC_SiGamma)
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

//================================================================================================
TEST(GPW_SCF, SolidCalculationMatchesTheSiAnchor)
{
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, ivec3_t(1,1,1));

    SCFParams par;
    par.NMaxIter=60; par.MinΔρ=1e-3; par.MinΔE=1e-6;
    par.MinΔFD=1e30; par.MinVirial=1e30; par.MinFD=1e30; par.StartingRelaxRo=0.3; par.MergeTol=1e-4;

    qchem::SolidCalculation calc(lat, MakeBasisSR(cell),
                                 {.Nelec=8, .species={{"Si",4}}, .densityEcut=20.0}, par);

    // N1/T1: the answers are reachable only through the PROOF, so a non-converged run cannot serve them.
    auto r = calc.Result();
    ASSERT_TRUE(r) << "SCF did not converge: " << (r ? std::string() : r.Error().details);
    EXPECT_NEAR(r->TotalCharge(), 8.0, 1e-6);
    EXPECT_NEAR(r->Energy(), -7.11506, 2e-3)        // the SiliconGammaConverges anchor, same tolerance
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
TEST(GPW_SCF, CrossRunFirstRunAnomalyProbe)
{
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, ivec3_t(1,1,1));

    SCFParams par;
    par.NMaxIter=60; par.MinΔρ=1e-3; par.MinΔE=1e-6;
    par.MinΔFD=1e30; par.MinVirial=1e30; par.MinFD=1e30; par.StartingRelaxRo=0.3; par.MergeTol=1e-4;

    double E[3];
    for (int r=0;r<3;++r)
    {
        std::cout<<"[xrun] ================= RUN "<<r<<" ================="<<std::endl;
        qchem::SolidCalculation calc(lat, MakeBasisSR(cell),
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
//  bit-for-bit, gated by GPW_SCF.CrossRunFirstRunAnomalyProbe above, so the arms run cold.)
static void ExpectRealComplexTwins(const Lattice_3D& lat, const UnitCell& cell, const char* what)
{
    const qchem::SolidCalcOptions optOn {.Nelec=8, .species={{"Si",4}}, .densityEcut=20.0};
    const qchem::SolidCalcOptions optOff{.Nelec=8, .species={{"Si",4}}, .densityEcut=20.0, .forceComplex=true};

    SCFParams one;                                   // the one-iteration recipe (gates 0 = strict-< never trips)
    one.NMaxIter=1; one.MinΔρ=0.0; one.MinΔE=0.0;
    one.MinΔFD=1e30; one.MinVirial=1e30; one.MinFD=1e30; one.StartingRelaxRo=0.3; one.MergeTol=1e-4;

    {   // (1) the arithmetic arm: exactly one iteration each
        qchem::SolidCalculation on(lat, MakeBasisSR(cell), optOn, one), off(lat, MakeBasisSR(cell), optOff, one);
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
        qchem::SolidCalculation on(lat, MakeBasisSR(cell), optOn, par), off(lat, MakeBasisSR(cell), optOff, par);
        auto ron=on.Result(), roff=off.Result();
        ASSERT_TRUE(ron)  << what << ": " << (ron  ? std::string() : ron .Error().details);
        ASSERT_TRUE(roff) << what << ": " << (roff ? std::string() : roff.Error().details);
        EXPECT_NEAR(ron->Energy(), roff->Energy(), 2e-5) << what;    // MinΔE=1e-6 resolution (measured 4.8e-6)
        const rvec3_t pts[]={ cell.ToCartesian(rvec3_t(0.3,0.4,0.7)),
                              cell.ToCartesian(rvec3_t(0.25,0.25,0.25)),
                              cell.ToCartesian(rvec3_t(0.1,0.9,0.2)) };
        for (const auto& r : pts)                                     // MinΔρ=1e-3 resolution (measured ≤9.1e-6)
            EXPECT_NEAR(ron->Density()(r), roff->Density()(r), 1e-4) << what;
    }
}

TEST(GPW_SCF, RealTRIMBlocksMatchComplex_SiGamma)
{
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, ivec3_t(1,1,1));   // Γ-only: EVERY block is TRIM, so the flip makes the whole run real
    ExpectRealComplexTwins(lat, cell, "Si Gamma real-vs-complex");
}

//  STEP 4's MIXED-MESH acceptance, at the smallest genuinely mixed mesh: N=(3,1,1) has ONE TRIM point
//  (Γ -- 2·1/3 is no reciprocal-lattice vector) beside two complex blocks (k=±1/3), so ONE composite
//  carries both child scalars through the full SCF -- the case Γ-only cannot reach.  (A 3×3×3 run is
//  the same code path 27 blocks wide; this keeps the gate's wall-time at ~3 Γ runs.)
TEST(GPW_SCF, RealTRIMBlocksMatchComplex_SiMixedMesh)
{
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, ivec3_t(3,1,1));
    ExpectRealComplexTwins(lat, cell, "Si (3,1,1) mixed-mesh real-vs-complex");
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
TEST(GPW_SCF, RealTRIMBlocksWithMOMMatchComplex_SiMixedMesh)
{
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, ivec3_t(3,1,1));

    SCFParams par;
    par.NMaxIter=60; par.MinΔρ=1e-3; par.MinΔE=1e-6;
    par.MinΔFD=1e30; par.MinVirial=1e30; par.MinFD=1e30; par.StartingRelaxRo=0.3; par.MergeTol=1e-4;
    par.UseMOM=true; par.MOMStartIter=2;              // armed EARLY: the reference must be live in-run

    qchem::SolidCalculation on (lat, MakeBasisSR(cell),
                                {.Nelec=8, .species={{"Si",4}}, .densityEcut=20.0}, par);
    qchem::SolidCalculation off(lat, MakeBasisSR(cell),
                                {.Nelec=8, .species={{"Si",4}}, .densityEcut=20.0, .forceComplex=true}, par);
    auto ron=on.Result(), roff=off.Result();
    ASSERT_TRUE(ron)  << "MOM on a mixed real/complex mesh must converge (R2.21): "
                      << (ron ? std::string() : ron.Error().details);
    ASSERT_TRUE(roff) << (roff ? std::string() : roff.Error().details);
    EXPECT_NEAR(ron->Energy(), roff->Energy(), 1e-9)
        << "a real TRIM block's MOM reference must behave exactly like its complex twin's";
    EXPECT_NEAR(ron->TotalCharge(), roff->TotalCharge(), 1e-10);
}


// The SHARP-FIELD leg of the gate (the plan names DISABLED_NaFRocksaltGamma as the stress case: the F-
// anion makes sharp peaks in rho and V_xc, and its diffuse basis is what the Becke grid exists for).
// DISABLED like the parent NaF anchor -- it is a long run (the full NaF convergence recipe at a
// matrix-grade densityEcut=160 reference, plus two Becke term evaluations); run it by hand with
// --gtest_also_run_disabled_tests when touching the XC quadrature.
TEST(GPW_SCF, BeckeXCMatchesUniformXC_NaFSR2)   // RE-ENABLED 2026-09-15: 24 s, Becke internally converged on the sharp-F system
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

// Mn PSEUDO-ATOM IN A BOX through the CRYSTAL (GPW) path vs the molecular facade -- the d-channel sibling
// of SiPseudoAtomInBoxMatchesFinite, and the cheap localiser for MnO's ~356 Ha over-binding.  Both
// real-space KB routes are now oracle-matched on this very PP (atomic radial -14.230 unpolarised /
// molecular Cartesian -14.668 polarised, vs CP2K ATOM -14.243986 restricted), so if the GPW path also
// lands near the facade the crystal KB is exonerated and the MnO defect lives elsewhere (multi-species,
// O, Ewald/alignment); if it does not, this is the l=2 crystal defect, isolated to ONE atom and one hour
// instead of the 4-atom magnetic cell.
// MEASURED 2026-08-06 -- THE ANSWER IS BASIS CONDITIONING, not the KB:
//   * CARTESIAN d carries the s contaminant (x^2+y^2+z^2), so our Mn window (7 s + 8 d shells x 6
//     Cartesian components = 55 functions) is RANK-DEFICIENT before any physics: lambdaMin 1.15e-07,
//     cond 8.2e7, and the GPW vet ABORTS.  The molecular facade only survives it by dropping modes --
//     its own log says "[ortho] dropped 5 near-null overlap mode(s) of 55".  A near-null direction that
//     the SCF can occupy is the classic variational-collapse mechanism, and MnO's 4-atom cell (154
//     functions, cond 7e8) sits in exactly that regime -- which fits -417 Ha vs the -61.47 oracle far
//     better than the 3e-2 analytic-vs-mesh KB discrepancy does.
//   * SPHERICAL d (5 pure components, no contaminant) is the natural cure but is NOT AVAILABLE on the
//     GPW path: it throws "the orbital basis is not a molecular Gaussian basis (no Gaussian::LatticeSum1E)"
//     -- the spherical lineage does not implement the lattice-sum face (cf. the parked S3b spherical work).
//   => CURED 2026-08-06 (user's insight): keep the d set and drop the s window to TWO functions.  The
//      contaminants already span the mid/tight s space, so only the DIFFUSE 4s tail (0.10) and one tight
//      s (24) are needed -- 2s+8d gives lambdaMin 3.0e-03 / cond 2.1e3 (from 1.15e-07 / 8.2e7) at a cost
//      of just 2 mHa (facade -14.6661 vs the rank-deficient 7s+8d's -14.6681).  Trimming the d count
//      instead "fixes" conditioning but costs 0.55 Ha (7s+4d -> -14.11): the d set is the physics, the
//      s window was the redundancy.  NB CP2K solves this same shell list SPHERICALLY (its log: 55
//      Cartesian vs 47 spherical functions) and never sees the contaminant -- apples-to-oranges when
//      comparing its oracles.  GPW_MN_SPHERICAL=1 re-runs the (throwing) spherical arm.
// THIS IS NOW A GATE: the first OCCUPIED-d species validated end to end through the crystal path.
TEST(GPW_SCF, MnAtomInBoxDChannel)
{
    // GPW_MN_SPHERICAL=1: the SPHERICAL arm (doc/SphericalLatticePlan.md I1) -- the facade reference then
    // runs the NATIVE spherical family (same span as the view), so the box-vs-facade A/B stays span-matched.
    const bool spherical=(bool)std::getenv("GPW_MN_SPHERICAL");
    const bool sphBasis =(bool)std::getenv("GPW_BASIS_SPH");     // the restored-s file (I3); spherical arms only
    Molecule mnmol; mnmol.Insert(new Atom(25, 0.0, {0,0,0}));
    Calculation cRef(mnmol, {.basis=sphBasis?"valence_lowq_sph":"valence_lowq_sr",
                             .multiplicity=6, .pseudopotential=true, .ppValence=7,
                             .angular = spherical ? Angular::Spherical : Angular::Cartesian});
    const double Eref=cRef.Energy();
    std::cout << "[Mn finite] "<<(sphBasis?"valence_lowq_sph":"valence_lowq_sr")<<" LSDA sextet (q7, "
              <<(spherical?"SPHERICAL":"CARTESIAN")
              <<")="<<Eref<<"   (CP2K ATOM UKS sextet oracle -14.674425)"<<std::endl;

    const Material box=qchem::Materials::Get("Mn_box16");
    const Lattice_3D lat=LatticeOf(box);
    SolidCalcOptions o=MnBoxOptions(box, "Mn atom-in-box sextet");   // S=5/2 Hund: nUp=6, nDown=1
    SCFParams par=Gates(40, 1e-5, 1e30); par.SmearingkT=5e-3;
    par.Verbose=(bool)std::getenv("GPW_MNO_VERBOSE");
    // CARTESIAN d carries the s CONTAMINANT (x^2+y^2+z^2), so 8 d shells duplicate the 7-function s space
    // -- measured lambdaMin 1.15e-07 / cond 8.2e7 on this one-atom box, i.e. the basis is rank-deficient
    // BEFORE any physics runs.  SPHERICAL d (5 pure components) removes the contaminant; GPW_MN_SPHERICAL=1
    // selects it for the A/B -- via the SPHERICAL LATTICE VIEW (doc/SphericalLatticePlan.md I1): the
    // NATIVE Angular::Spherical family has no LatticeSum1E capability (the historical blocker -- feeding
    // it here died on the GPW cross-cast), so the view over the Cartesian engine is the working door.
    std::shared_ptr<const Real_BS> mnbasis(
        BasisSet::Gaussian::Factory(sphBasis?BasisSetData::VALENCE_LOWQ_SPH:BasisSetData::VALENCE_LOWQ_SR,
                                    box.cell.get(), BasisSet::Gaussian::Engine::MnD,
                                    BasisSet::Gaussian::Angular::Cartesian));
    if (spherical) mnbasis=BasisSet::Gaussian::PG_Spherical::MakeSphericalLatticeView(mnbasis);
    std::cout << "[Mn in-box] angular=" << (spherical?"SPHERICAL":"CARTESIAN") << std::endl;
    GpwReport report("Mn "+o.label, par.Verbose);
    qchem::SolidCalculation calc(lat, mnbasis, o, par);
    const double E=calc.LastIterateTerms().GetTotalEnergy();   // as before: the energy is pinned, convergence is not asserted
    std::cout << "[Mn in-box] GPW="<<E<<"  facade="<<Eref<<"  diff="<<(E-Eref)<<std::endl;
    EXPECT_NEAR(calc.LastIterateCharge(), 7.0, 1e-6);
    EXPECT_NEAR(E, Eref, 3e-2) << "GPW d-channel vs the molecular facade (measured 12 mHa)";
    if (!spherical)
        EXPECT_NEAR(E, -14.6380, 1e-3);                      // did-E-move anchor (2s+7d, the SR-trimmed cell basis)
}

// (PolarizedRunKeepsItsSpin -- the Mn sextet asking for Kerker, 217 s, 27% of the suite -- was DELETED
//  2026-09-15 (doc/TestSuitePlan.md phase 1).  Its claim is a MIXER property: a polarized seed must get a
//  per-channel ρ̃ mixer, never a single-map one.  It is now pinned with no SCF at all by
//  src/ChargeDensity/tests/KerkerMix.C -- PolarizedSeedComposesPerChannel, PolarizedStepMovesEachChannelAtAlpha
//  and the QCHEM_SPINBLIND_KERKER negative control -- in 41 ms.  What else it asserted is covered by
//  MnAtomInBoxDChannel above (same cell, the polarized energy) and ImposedShubnikovHoldsAFMThroughSCF_Mn2Box
//  below (order SUSTAINED through an SCF).)

// ============================ MnO rocksalt AFM-II (SymmetryUpgradePlan §7 step 7) ============================
// The campaign RUN (RunMnO: the rhombohedral AFM-II cell, the recipe, the anneal schedule, the FM arm and the
// ordering comparison) is CLIapps/gpwprobe's `mno` sub-command since 2026-09-15 (doc/TestSuitePlan.md §8).  What
// stays here are the SEED-LEVEL gates: no SCF, cheap, and the mirror they pin is exact.

// THE RAW SEED of the MnO AFM-II cell: are the two Mn sublattices actually equal and opposite?
//
// One magnetic species on one Wyckoff site means the two Mn are related by the (1/2,1/2,1/2) translation,
// so ANY valid magnetic solution -- seed included -- must have m(Mn1) = -m(Mn2).  Nothing checked that.
// PlaneWaveDFT.PolarizedSeedAFMStaggering covers a DIFFERENT cell (simple cubic, 2 Mn, no O, NEUTRAL
// targets) and asserts G-space quantities plus GetTotalSpin()==0; m_stag = 1/2(m1-m2), the campaign's order
// parameter, is BLIND to the imbalance -- it reads +0.366 for (+0.37,-0.37) and for (0,-0.73) alike.
// Measured 2026-08-11: the density one Fock build downstream of this seed has m1 = -0.0001, m2 = -0.73.
// This test asks whether the seed ITSELF is lopsided, i.e. whether the defect is upstream of the SCF.
TEST(GPW_SCF, MnOSeedSublatticesAreEqualAndOpposite)
{
    const double a=8.40;
    const Material mno=qchem::Materials::Get("MnO_AFM2");   // the rhombohedral AFM-II cell: Mn +m at 0 (the CORNER atom), Mn -m at 1/2, O at 1/4, 3/4
    const std::shared_ptr<UnitCell> cellp=mno.cell;
    UnitCell& cell=*cellp;
    Lattice_3D lat(cell, ivec3_t(1,1,1));

    // The seed the run actually uses: PolarizedSeedCD over this cell, with the IonicSAD targets
    // (Mn2+ => 5 of the q7 valence, O2- => 8 of the q6).  Built directly -- no SCF, no Hamiltonian.
    qchem::BasisSet::PlaneWave::PlaneWave_IBS pw(lat.Reciprocal(), ivec3_t(1,1,1), ivec3_t(0,0,0), 4.0);
    std::shared_ptr<const BasisSet::cFIT_CD_ABS> fb(pw.CreateCDFitBasisSet(&cell, qcMesh::MeshParams{}));
    const std::map<size_t,int> ionic{{25,5},{8,8}};
    qchem::ChargeDensity::PolarizedSeedCD seedCD(fb, &cell, "LDA", ionic);
    const auto* up=seedCD.GetChannel(Spin::Up);
    const auto* dn=seedCD.GetChannel(Spin::Down);

    const rvec3_t off(0.7,0,0), rMn1(0,0,0), rMn2(a,a,a);   // A*(1/2,1/2,1/2) = a(1,1,1); 0.7 bohr = the d peak
    const double m1=(*up)(rMn1+off)-(*dn)(rMn1+off);
    const double m2=(*up)(rMn2+off)-(*dn)(rMn2+off);
    std::cout << "[seed probe] m1=" << m1 << " m2=" << m2 << " m1+m2=" << m1+m2
              << "  N_up=" << up->GetTotalCharge() << " N_dn=" << dn->GetTotalCharge() << std::endl;

    EXPECT_NEAR(up->GetTotalCharge(), dn->GetTotalCharge(), 1e-8) << "the AFM seed carries no NET moment";
    EXPECT_GT(std::abs(m1), 0.1) << "the + sublattice must actually be polarized";
    EXPECT_GT(std::abs(m2), 0.1) << "the - sublattice must actually be polarized";
    EXPECT_NEAR(m1+m2, 0.0, 0.02*std::abs(m1))
        << "the two sublattices are related by a lattice translation: m1 must equal -m2";

    // THE MIRROR RELATION EVERYWHERE, not just at the two nuclei.  The AFM seed satisfies
    // rho_up(r) = rho_dn(r + t) with t = A*(1/2,1/2,1/2) = (a,a,a) at EVERY r, by construction.  Probed
    // here at points that straddle the CELL BOUNDARY, because that is where the Becke XC mesh wraps its
    // points into the home cell (kpt = r - A*n0) and where a real-space evaluation built from per-atom
    // recentred radials WITHOUT periodic images would pick up the wrong atom -- which, the two Mn carrying
    // OPPOSITE spins, swaps the channels rather than merely losing density.  That failure is invisible to
    // the grid-charge check (rho_up+rho_dn is untouched; only rho_up-rho_dn flips) and it is exactly the
    // observed symptom: the corner atom's moment dies while the mesh integrates the total to 1e-5.
    struct P { const char* what; rvec3_t r; };
    const std::vector<P> probes = {
        {"inside, +x of Mn1",        rvec3_t( 0.7, 0.0, 0.0)},
        {"OUTSIDE the cell, -x",     rvec3_t(-0.7, 0.0, 0.0)},   // physically beside Mn1; wraps
        {"OUTSIDE the cell, -xyz",   rvec3_t(-0.5,-0.5,-0.5)},
        {"far tail, -2 bohr",        rvec3_t(-2.0, 0.0, 0.0)},
    };
    for (const auto& p : probes)
    {
        const rvec3_t rt = p.r + rvec3_t(a,a,a);                  // the mirror partner
        const double u=(*up)(p.r), d=(*dn)(rt);
        std::cout << "[mirror] " << p.what << ": rho_up(r)=" << u << "  rho_dn(r+t)=" << d
                  << "  diff=" << u-d << std::endl;
        EXPECT_NEAR(u, d, 1e-6*std::max(1.0,std::abs(u)))
            << "AFM mirror broken at " << p.what << ": the seed is not properly periodic there";
    }

    // THE BATCH OVERLOAD vs THE SINGLE-POINT ONE.  Everything above used operator()(rvec3_t).  The XC mesh
    // samples a MATRIX-FREE seed through the BATCHED operator()(rvec3vec_t) instead -- SinglesDensitySampler::RhoPol
    // takes its cSpinResolved_CD branch for exactly this density -- so the batch path is what the first Fock
    // build actually sees, and nothing has ever checked the two agree.  They must, pointwise.
    rvec3vec_t batch(2*probes.size());
    for (size_t i=0;i<probes.size();++i) { batch[i]=probes[i].r; batch[probes.size()+i]=probes[i].r+rvec3_t(a,a,a); }
    const rvec_t bu=(*up)(batch), bd=(*dn)(batch);
    ASSERT_EQ(bu.size(), batch.size());
    for (size_t i=0;i<batch.size();++i)
    {
        const double su=(*up)(batch[i]), sd=(*dn)(batch[i]);
        std::cout << "[batch] r=(" << batch[i].x << "," << batch[i].y << "," << batch[i].z << ")"
                  << " up: batch=" << bu[i] << " single=" << su << " d=" << bu[i]-su
                  << " | dn: batch=" << bd[i] << " single=" << sd << " d=" << bd[i]-sd << std::endl;
        EXPECT_NEAR(bu[i], su, 1e-6*std::max(1.0,std::abs(su))) << "UP batch != single at probe " << i;
        EXPECT_NEAR(bd[i], sd, 1e-6*std::max(1.0,std::abs(sd))) << "DN batch != single at probe " << i;
    }
}

// THE v_xc SUBLATTICE MIRROR ON THE XC MESH -- the SymmetryUpgradePlan "NEXT ACTION" probe (2026-08-11).
// By elimination (seed, Becke weights, Kinetic/Vloc/Vnl, Phi tables all exonerated) the first-Fock-build
// mirror break must live in v_xc -- yet a pointwise LSDA functional "cannot" be site-dependent.  This probe
// resolves the contradiction by testing what the Fock build ACTUALLY consumes: the channel rasters
// SinglesDensitySampler::RhoPol hands the spin-native XC term, at the mesh's own points.  The mesh stores its
// points WRAPPED into the home cell (kpt = r - A*n0, MakePeriodicBeckeMesh) -- so a valid seed must satisfy
// rho_up(p_g) = rho_dn(p_g + t) with t = A*(1/2,1/2,1/2) AT EVERY STORED POINT, and (v_xc being pointwise
// in the channel pair) v_xc^up(p_g) = v_xc^dn(p_g + t).  The two Mn blocks' grids are exact t-translates of
// each other (same radial x angular template), so every point's mirror partner is itself a mesh point --
// found here by hashing wrapped fractional coordinates, no interpolation anywhere.
//
// WHY THE SEED GATE ABOVE COULD PASS WHILE THIS FAILS: its probe points sit within ~2 bohr of a HOME atom,
// where SeedCD's real-space rho(r) = Sum_atoms rho_atom(|r-R|) -- a sum with NO LATTICE IMAGES -- is
// dominated by an atom it actually contains.  A WRAPPED mesh point near a cell face reads its density from
// an IMAGE of an atom (for the CORNER Mn, 7 of the 8 octants of its density hump belong to images), which
// an image-less sum simply does not have.  That is site-specific (the centre Mn2 has no near-shell wrapped
// points), a rigid translation changes it (MNO_SHIFT), the uniform raster reproduces it (its corner-region
// points read image density too), and under the AFM staggering the missing hump is the MAJORITY channel of
// exactly one sublattice -- every recorded symptom.
TEST(GPW_SCF, MnOSeedVxcMirrorOnBeckeMesh)
{
    const double a=8.40;
    const Material mno=qchem::Materials::Get("MnO_AFM2");   // the rhombohedral AFM-II cell: Mn +m at 0 (the CORNER atom), Mn -m at 1/2, O at 1/4, 3/4
    const std::shared_ptr<UnitCell> cellp=mno.cell;
    UnitCell& cell=*cellp;
    Lattice_3D lat(cell, ivec3_t(1,1,1));

    // The seed the run uses (identical to MnOSeedSublatticesAreEqualAndOpposite).
    qchem::BasisSet::PlaneWave::PlaneWave_IBS pw(lat.Reciprocal(), ivec3_t(1,1,1), ivec3_t(0,0,0), 4.0);
    std::shared_ptr<const BasisSet::cFIT_CD_ABS> fb(pw.CreateCDFitBasisSet(&cell, qcMesh::MeshParams{}));
    const std::map<size_t,int> ionic{{25,5},{8,8}};
    qchem::ChargeDensity::PolarizedSeedCD seedCD(fb, &cell, "LDA", ionic);

    // The run's Becke recipe, shrunk (nR=20, GL-11) -- wrapping is generic, it needs no production density.
    qcMesh::MeshParams mp=qcMesh::BeckeXCParams(20, -1.0, 11);
    auto mesh=std::make_shared<const qcMesh::Mesh>(cell.CreateIntegrationMesh(mp));
    const rvec3vec_t& P=mesh->Points();
    const size_t N=P.size();
    ASSERT_GT(N, 0u);

    // The channel rasters EXACTLY as the Fock build gets them (RhoPol's cSpinResolved_CD seed branch).
    auto engine=SinglesEngineOver({mesh, {}});
    const rvec_t ru=engine->RhoPol(&seedCD, Spin::Up);
    const rvec_t rd=engine->RhoPol(&seedCD, Spin::Down);

    // Mirror-partner lookup: hash each point's wrapped fractional coords; partner(g) = the mesh index at
    // wrapped(frac(p_g) + (1/2,1/2,1/2)).  Quantized key + 27-neighbour probe rides out wrap roundoff.
    const double q=1e-9;
    auto wrapf=[](rvec3_t f){ f.x-=floor(f.x); f.y-=floor(f.y); f.z-=floor(f.z); return f; };
    std::map<std::tuple<long long,long long,long long>,size_t> at;
    std::vector<rvec3_t> F(N);
    for (size_t g=0; g<N; g++)
    {
        F[g]=wrapf(cell.ToFractional(P[g]));
        at[{llround(F[g].x/q),llround(F[g].y/q),llround(F[g].z/q)}]=g;
    }
    auto partner=[&](size_t g)->long
    {
        const rvec3_t fm=wrapf(F[g]+rvec3_t(0.5,0.5,0.5));
        const long long kx=llround(fm.x/q), ky=llround(fm.y/q), kz=llround(fm.z/q);
        for (long long dx=-1; dx<=1; dx++) for (long long dy=-1; dy<=1; dy++) for (long long dz=-1; dz<=1; dz++)
            if (auto it=at.find({kx+dx,ky+dy,kz+dz}); it!=at.end())
            {
                rvec3_t d=wrapf(F[it->second]-fm+rvec3_t(0.5,0.5,0.5))-rvec3_t(0.5,0.5,0.5);  // min-image
                if (norm(d)<5e-9) return long(it->second);
            }
        return -1;
    };

    // v_xc per point from the channel pair -- the same functionals the spin-native XC term applies, through
    // the same two-channel face.
    qchem::Hamiltonian::SlaterExchange ex(2.0/3.0);
    qchem::Hamiltonian::VWN_Correlation vc;
    auto vxc=[&](double u, double d, const Spin& s)->double
    {
        if (u+d<=1e-12) return 0.0;                       // VWN's r_s/log guard; symmetric, mirror-safe
        return ex.GetVxc(u, d, s) + vc.GetVxc(u, d, s);
    };

    // Sweep: rho and v_xc mirror defects across the whole raster; localize the worst offenders.
    size_t nOrphan=0, nBad=0;
    double maxRho=0, maxV=0, maxOrphanWRho=0;
    size_t argRho=0; long argRhoJ=-1;
    const rvec_t& W=mesh->Weights();
    for (size_t g=0; g<N; g++)
    {
        const long j=partner(g);
        // An ORPHAN is legal but must be an eps-tail point: the free Becke builder's keep decisions are
        // bit-different between translated blocks, so an eps-BORDERLINE point can be kept on one atom and
        // dropped on its partner (the same mechanism the imposed builder's orbit-consistency pass drops).
        // What the quadrature sees of it is w*rho -- bound THAT, not the count.
        if (j<0) { nOrphan++; maxOrphanWRho=std::max(maxOrphanWRho, std::abs(W[g])*std::max(ru[g],rd[g])); continue; }
        const double dRho=std::max(std::abs(ru[g]-rd[j]), std::abs(rd[g]-ru[j]));
        const double dV  =std::max(std::abs(vxc(ru[g],rd[g],Spin::Up  )-vxc(ru[j],rd[j],Spin::Down)),
                                   std::abs(vxc(ru[g],rd[g],Spin::Down)-vxc(ru[j],rd[j],Spin::Up  )));
        if (dRho>1e-6) nBad++;
        if (dRho>maxRho) { maxRho=dRho; argRho=g; argRhoJ=j; }
        if (dV  >maxV) maxV=dV;
    }
    std::cout << "[vxc mirror] N=" << N << " orphans=" << nOrphan << " (max w*rho " << maxOrphanWRho
              << ") bad(rho>1e-6)=" << nBad
              << "  max|rho_up(r)-rho_dn(r+t)|=" << maxRho
              << "  max|vxc_up(r)-vxc_dn(r+t)|=" << maxV << std::endl;

    // Localize the worst point: whose density is it -- a HOME atom's, or an IMAGE's the seed cannot see?
    if (argRhoJ>=0)
    {
        std::vector<rvec3_t> R; std::vector<std::string> nm={"Mn1","Mn2","O1","O2"};
        for (auto atom : cell) R.push_back(atom->itsR);
        auto nearest=[&](const rvec3_t& r)
        {
            double best=1e300; std::string who;
            for (size_t ia=0; ia<R.size(); ia++)
                for (int i=-1;i<=1;i++) for (int jj=-1;jj<=1;jj++) for (int k=-1;k<=1;k++)
                {
                    const double d=norm(r-(R[ia]+cell.ToCartesian(rvec3_t(i,jj,k))));
                    if (d<best) { best=d; who=nm[ia]+((i||jj||k)?" IMAGE":""); }
                }
            return std::make_pair(best,who);
        };
        const auto [dg,wg]=nearest(P[argRho]);
        const auto [dj,wj]=nearest(P[size_t(argRhoJ)]);
        std::cout << "[vxc mirror] worst point r=("<<P[argRho].x<<","<<P[argRho].y<<","<<P[argRho].z
                  << ") nearest "<<wg<<" d="<<dg<<"  rho_up(r)="<<ru[argRho]<<" rho_dn(r)="<<rd[argRho]<<"\n"
                  << "[vxc mirror] partner     r=("<<P[size_t(argRhoJ)].x<<","<<P[size_t(argRhoJ)].y<<","
                  << P[size_t(argRhoJ)].z<<") nearest "<<wj<<" d="<<dj
                  << "  rho_up(r+t)="<<ru[size_t(argRhoJ)]<<" rho_dn(r+t)="<<rd[size_t(argRhoJ)]<<std::endl;
    }

    EXPECT_LT(maxOrphanWRho, 1e-8) << "a partnerless mesh point must be an eps-tail point (weight*rho "
                                      "below the Becke builder's eps-converged-series contract)";
    EXPECT_LT(maxRho, 1e-8) << "the seed's channel rasters break the sublattice mirror ON THE XC MESH -- "
                               "this is the site-dependent v_xc defect (SymmetryUpgradePlan WHERE-WE-LEFT-OFF)";
    EXPECT_LT(maxV, 1e-6) << "v_xc^up(r) != v_xc^dn(r+t) on the mesh: the first Fock build is fed a "
                             "mirror-broken exchange-correlation potential";
}

// SHUBNIKOV S3, end to end at the seed level (doc/SymmetryUpgradePlan.md §7 step 7): a MAGNETICALLY
// imposed MnO basis must star-average the channel pair under the Shubnikov group -- keeping the AFM
// staggering EXACTLY mirror-symmetric -- where the grey average would erase it.  Chain under test:
// MagneticDecoration (the seed's own species rule) -> GPWParams::siteSpins -> the factory's Shubnikov
// resolution -> CreateXCQuadrature (site-adapted invariant mesh + fold + sigma tags + flip-fixed zero
// flags) -> SinglesDensitySampler::RhoPol's (rho,m) split.
TEST(GPW_SCF, MnOImposedShubnikovKeepsTheSeedStaggering)
{
    namespace L3=BasisSet::Lattice;
    using namespace qchem::ChargeDensity;
    const double a=8.40;
    const Material mno=qchem::Materials::Get("MnO_AFM2");   // the rhombohedral AFM-II cell: Mn +m at 0 (the CORNER atom), Mn -m at 1/2, O at 1/4, 3/4
    const std::shared_ptr<UnitCell> cellp=mno.cell;
    UnitCell& cell=*cellp;
    Lattice_3D lat(cell, ivec3_t(1,1,1));

    // The decoration, by the seed's own resolution: Mn2+ (d5 pair) = +/-1, O2- (closed shell) = 0.
    auto st=lat.GetStructure();
    std::vector<int> spins=MagneticDecoration(st.get(), "LDA", IonicSADTargets(st.get(), "LDA"));
    ASSERT_EQ(spins.size(), 4u);
    EXPECT_EQ(spins[0], +1); EXPECT_EQ(spins[1], -1);
    EXPECT_EQ(spins[2],  0); EXPECT_EQ(spins[3],  0);

    // The magnetically IMPOSED GPW basis (coarse everything: this gate tests symmetry, not accuracy).
    std::unique_ptr<Complex_BS> bs(L3::GPWFactory(lat, MakeBasisLowQ(cell, BasisSetData::VALENCE_LOWQ_SR),
        L3::GPWParams{.densityEcut=8.0, .imposeSymmetry=true, .siteSpins=spins}));
    BasisSet::FitQuadrature q = bs->CreateXCQuadrature(st.get(), qcMesh::BeckeXCParams(20, -1.0, 11));
    ASSERT_EQ(q.NumSpinOps(), 24u) << "the Shubnikov group of the AFM-II cell has 24 ops (12+12)";
    size_t nFlip=0; for (auto s : q.GetSpinOps()) if (s==Symmetry::SpinAction::Flip) nFlip++;
    EXPECT_EQ(nFlip, 12u);
    ASSERT_EQ(q.NumFlipFixed(), q.GetMesh()->size());   // (FoldedMesh's ctor now checks this too)

    // The seed (the same PolarizedSeedCD the run uses), through the engine's channel-pair projector.
    qchem::BasisSet::PlaneWave::PlaneWave_IBS pw(lat.Reciprocal(), ivec3_t(1,1,1), ivec3_t(0,0,0), 4.0);
    std::shared_ptr<const BasisSet::cFIT_CD_ABS> fb(pw.CreateCDFitBasisSet(&cell, qcMesh::MeshParams{}));
    const std::map<size_t,int> ionic{{25,5},{8,8}};
    PolarizedSeedCD seedCD(fb, &cell, "LDA", ionic);

    auto engine=SinglesEngineOver(q);
    const rvec_t up=engine->RhoPol(&seedCD, Spin::Up);
    const rvec_t dn=engine->RhoPol(&seedCD, Spin::Down);

    // (a) The staggering SURVIVES the imposed projector: the magnetization raster keeps its full scale.
    double maxM=0; for (size_t g=0; g<up.size(); g++) maxM=std::max(maxM, std::abs(up[g]-dn[g]));
    EXPECT_GT(maxM, 0.1) << "the SHUBNIKOV average must PRESERVE the AFM staggering";

    // (b) ...and is EXACTLY mirror-antisymmetric on the raster: m must vanish at every flip-fixed point,
    // and the weighted net moment must be zero to machine precision (the signed projector's guarantees).
    const rvec_t& W=q.GetMesh()->Weights();
    double net=0, worstFixed=0;
    for (size_t g=0; g<up.size(); g++)
    {
        net += W[g]*(up[g]-dn[g]);
        if (q.IsFlipFixed(g)) worstFixed=std::max(worstFixed, std::abs(up[g]-dn[g]));
    }
    // The projector pairs every orbit's members +/- exactly; the residual is the ULP-level asymmetry of
    // the site-adapted builder's partner WEIGHTS (op-image copies, equal only to roundoff -- measured
    // 1.7e-9 over 6000 points), not a projector defect.
    EXPECT_LT(std::abs(net), 1e-7) << "an imposed AFM pair carries ZERO net moment by construction";
    EXPECT_LT(worstFixed, 1e-12) << "m must vanish exactly at the flip-fixed mesh points";

    // (c) The GREY control: the SAME mesh and fold with the sigma tags withheld = the historical
    // per-channel spatial average, which maps +m sites onto -m sites and ERASES the order.  This is the
    // unit-level half of the S4 negative control ("imposing the grey group kills m_stag").
    auto grey=SinglesEngineOver(qcMesh::FoldedMesh(q.GetMesh(), q.GetFold()));
    const rvec_t gup=grey->RhoPol(&seedCD, Spin::Up);
    const rvec_t gdn=grey->RhoPol(&seedCD, Spin::Down);
    double maxGrey=0; for (size_t g=0; g<gup.size(); g++) maxGrey=std::max(maxGrey, std::abs(gup[g]-gdn[g]));
    EXPECT_LT(maxGrey, 1e-10) << "the grey average must ERASE the staggering -- if it does not, the "
                                 "Shubnikov machinery is not actually load-bearing";

    // (d) The TOTAL density is identical through both engines (the even channel is sigma-blind).
    double dTot=0; for (size_t g=0; g<up.size(); g++) dTot=std::max(dTot, std::abs((up[g]+dn[g])-(gup[g]+gdn[g])));
    EXPECT_LT(dTot, 1e-10) << "sigma must not touch the total density";
}

// SHUBNIKOV S4, THROUGH SCF (doc/SymmetryUpgradePlan.md §7 step 7): the imposed magnetic star-average
// confines the magnetization to the STAGGERED sector by construction -- every iterate's m is projected
// onto the Shubnikov-symmetric cone, so the AFM order cannot leak into a net moment and the sublattice
// mirror holds EXACTLY at every iteration, converged or not.  Fixture: the cheapest genuinely staggered
// crystal -- two neutral Mn (the d5s2 Hund pair) in a cubic box, CsCl/B2 arrangement, AFM flip on the
// second.  Bounded iterations (this gate tests SYMMETRY through the live SCF loop, not convergence).
// The through-SCF GREY negative control lives on MnO (MNO_IMPOSE=2): THIS cell's detected grey group is
// all site-preserving (every cubic W admits tau=0), so grey imposition here would be a vacuous control.
TEST(GPW_SCF, ImposedShubnikovHoldsAFMThroughSCF_Mn2Box)
{
    const Material box=qchem::Materials::Get("Mn2_box7");     // Mn +m at 0, Mn -m at 1/2 (the AFM flip), a=7
    const double a=7.0;
    const Lattice_3D lat=LatticeOf(box);
    SolidCalcOptions o=OptionsFor(box, "Mn2 B2 AFM imposed");
    o.multiplicity=1;
    o.seed=qchem::ChargeDensity::SeedStrategy::SAD;      // neutral Mn: the library's d5s2 spin pair
    o.imposeSymmetry=true;                                // S3 resolves the SHUBNIKOV group from the flips
    o.ortho=qchem::CholeskyPivoted; o.orthoTol=1e-4;
    o.xcMesh=qcMesh::BeckeXCParams(20, -1.0, 11);         // coarse quadrature: symmetry, not accuracy
    o.xcMesh.cellKind=qcMesh::UnitCellKind::Becke;
    SCFParams par=Gates(6, 1e-9, 1e30);                   // BOUNDED: this gate tests symmetry through the loop, not convergence
    par.SmearingkT=5e-3;                                  // the open-d-manifold tie smoother
    par.StartingRelaxRo=0.45;
    par.Verbose=(bool)std::getenv("GPW_MN2_VERBOSE");
    GpwReport report("Mn "+o.label, par.Verbose);
    qchem::SolidCalculation calc(lat, MakeBasisLowQ(*box.cell, BasisSetData::VALENCE_LOWQ_SR), o, par);
    const GpwHandles h=LastIterateHandles(calc);
    ASSERT_TRUE(h.cd) << "the run must produce a final density";

    // The final density through the spin-resolved face: the order must be ALIVE and EXACTLY mirrored.
    const auto* pol=dynamic_cast<const qchem::ChargeDensity::cSpinResolved_CD*>(h.cd);
    ASSERT_NE(pol, nullptr);
    const auto* up=pol->GetChannel(Spin::Up);
    const auto* dn=pol->GetChannel(Spin::Down);
    ASSERT_TRUE(up && dn);
    const rvec3_t off(0.7,0,0), r1(0,0,0), r2(a/2,a/2,a/2);
    const double m1=(*up)(r1+off)-(*dn)(r1+off);
    const double m2=(*up)(r2+off)-(*dn)(r2+off);
    std::cout << "[Mn2 imposed] after "<<calc.IterationCount()<<" iterations: m1="<<m1<<" m2="<<m2
              << " m1+m2="<<m1+m2<<" Etot="<<calc.LastIterateTerms().GetTotalEnergy()<<std::endl;
    EXPECT_GT(std::abs(m1), 0.05) << "the AFM order must SURVIVE the imposed SCF loop";
    // What the projector holds EXACTLY mirrored is the rho the FOCK consumes (the engine's (rho,m)
    // star-average); the DENSITY MATRIX itself is deliberately NOT projected (the rho-projection
    // philosophy, T3.2 note), so the D-density probed here is DRIVEN toward the mirror at the
    // convergence rate -- measured |m1+m2| = 3.1e-3 beside lastΔρ = 5.6e-4 at the 6-iteration cap.
    // The gate therefore asserts the mirror at the "small beside the order" tier; the EXACT-mirror
    // guarantee on the projected rasters is the S3 seed-level gate above.
    EXPECT_NEAR(m1+m2, 0.0, 0.05*std::abs(m1))
        << "the imposed SCF must keep the sublattice mirror tight (driven by the projected Fock)";
    EXPECT_NEAR(pol->GetTotalSpin(), 0.0, 1e-8) << "an imposed AFM pair carries zero net moment";
}

// ★ T2's POSITIVE PATH (doc/OpenWork.md N1/T2, built 2026-08-26).  The OrderLost postcondition shipped
// with only its negative side exercised: every order-losing run measured on 2026-08-25 ALSO ran out of
// iterations, so NotConverged tripped first and the branch that MINTS an OrderLost had never once run.
// An unexercised guard is a guess, so this is the case that fires it.
//
// THE FIXTURE, and why each ingredient is load-bearing:
//   - Na2 in a box at its bond length.  Two neutral Na, so the SAD seed plants the library's spin pair on
//     each site (+-1 e), while the GROUND STATE is the closed-shell sigma^2 singlet -- m=0 is the honest
//     answer here, which is exactly the contradiction T2 exists to catch when it is asked for beside an
//     imposition.  It is also the cheapest such cell: two valence electrons.
//   - AFM FLIP on the second atom -- without it the decoration is ferromagnetic and the seed carries no
//     STAGGERING for the imposed Shubnikov group to preserve.
//   - imposeSymmetry -- T2 is gated on it (a FREE run finding m=0 is physics and must never be touched).
//   - A BECKE mesh, PINNED.  The postcondition reads the INTEGRATED site moment, and the integration
//     basins are the Becke site blocks; a uniform mesh has none and SiteMoments correctly returns empty,
//     so the check silently skips.  Auto would route this soft cell to the uniform mesh.
// The mirror-image gate is ImposedShubnikovHoldsAFMThroughSCF_Mn2Box above: same shape, real magnet,
// order SURVIVES.  The two together say the detector discriminates rather than always firing.
TEST(GPW_SCF, ImposedOrderLostIsAPostconditionFailure_Na2Box)
{
    // ⛔ THIS GATE CANNOT RUN WITHOUT THE BECKE MESH, and the header above already says why: the
    // postcondition reads the INTEGRATED SITE MOMENT, whose basins ARE the Becke site blocks.  Vetoing the
    // atom-centred mesh (QCHEM_BECKE_XC=0 / CP2K_COMPAT=1) leaves SiteMoments correctly empty, the run is
    // then a perfectly good non-magnetic answer, and OrderLost has nothing to fire on -- the diagnostics
    // say so in as many words ("order: not measurable (no atom-centred basins on this XC mesh)").  So the
    // honest outcome is SKIP: the subject is unmeasurable in that configuration, not broken by it.
    if (!qchem::theRunPolicy().BeckeXC())
        GTEST_SKIP() << "the OrderLost postcondition is measured on Becke site basins, which the run "
                        "policy has vetoed (QCHEM_BECKE_XC=0 / CP2K_COMPAT=1)";
    const double a=16.0, d=5.8;                       // ~Na2 bond length (au) in the Si gate's box
    UnitCell cell(a);
    cell.AddAtom(11, {0.5-0.5*d/a,0.5,0.5}, false);   // Na +m
    cell.AddAtom(11, {0.5+0.5*d/a,0.5,0.5}, true);    // Na -m -- the AFM flip the SAD seed plants
    Lattice_3D lat(cell, ivec3_t(1,1,1));

    SCFParams par;
    // alpha=0.5 and a generous cap are LOAD-BEARING, not decoration: this gate needs a run that
    // CONVERGES, because a run that merely hits the cap trips NotConverged first and leaves the
    // postcondition unexercised all over again.
    //
    // ⚠ THIS WAS A FLAKY TEST, and the 2026-08-27 sweep says why (user: "textbook example of a flaky
    // test ... much more flaky than average+3sigma").  It binds a DETERMINISTIC subject -- does
    // OrderLost get minted -- to a NON-deterministic precondition: whether a marginal SCF converges
    // inside the cap.  A 100-iteration cap left 34 iterations of margin, and a change four layers down
    // in the collocation kernel spent all of it.  The detector RULES are now unit-tested on synthetic
    // trajectories (src/Calculation/tests/RunDiagnostics.C), so this gate is no longer the only thing
    // standing between us and an untested postcondition; what remains here is the WIRING.
    //
    // MEASURED 2026-08-27, sweeping alpha x criterion x kernel (NA2_* knobs below), NMaxIter=400:
    //   alpha  0.3: 311 iters (walk) / 306 (contracted)      <- NOT "oscillates forever"; the earlier
    //   alpha  0.4: 300              / never (oscillates)       comment said that against a 100 cap
    //   alpha  0.5:  66              / 290                    <- 0.5+walk is a lucky sweet spot
    // and five of those six share ONE mode: E is dead at -0.3320448527 (dE~1e-15) by iteration ~200
    // while Δρ grinds down geometrically, every run stopping the instant it crosses MinΔρ.  So the long
    // tail is not oscillation, it is a slow LINEAR density mode.
    // ⇒ SMALLER ALPHA DOES NOT HELP (0.3 is slower than 0.5, and 0.4 is worse than either).
    // ⇒ NEITHER DOES DENSITY-PULAY, and its failure is the instructive one: NA2_PULAY=5 converges in 8
    //   iterations (walk) / 14 (contracted) -- but to a DIFFERENT state, E 4.4e-5 higher with the moment
    //   still alive at 0.031 e, so OrderLost never fires and the gate is defeated a third way.  The
    //   SCFParams note about history mixing being unstable far from the fixed point applies here.
    // ⇒ THE LONG RUN IS NOT SLOP: the fixture NEEDS it, because the MOMENT is a slower mode than Δρ and
    //   has to fall below 1% of its peak before the postcondition can fire at all.
    // So the cap is raised to 400 -- the minimal honest change, buying margin without altering what is
    // tested.  Cost: ~3 s on the walk, ~13 s on the contracted kernel.  Unwinding an AFM seed that the
    // answer does not want is simply harder than starting unpolarized: the same cell unpolarized
    // converges in 30.
    // NA2_ALPHA / NA2_NMAX / NA2_TRACE: the fixture's convergence recipe is INSTRUMENTAL (it exists so
    // the run converges, so NotConverged does not pre-empt the OrderLost this gate is named for), which
    // makes it exactly the thing an investigation needs to sweep.  Same idiom as MNO_ALPHA.
    auto envd2=[](const char* n, double d){ const char* v=std::getenv(n); return v ? std::atof(v) : d; };
    auto envi2=[](const char* n, int    d){ const char* v=std::getenv(n); return v ? std::atoi(v) : d; };
    par.NMaxIter=envi2("NA2_NMAX",400); par.MinΔρ=envd2("NA2_DRHO",1e-6); par.MinΔE=1e30;
    par.MinΔFD=1e30; par.MinVirial=1e30; par.MinFD=1e30;
    par.MergeTol=1e-4;
    par.PulayDepth=envi2("NA2_PULAY",0); par.PulayStart=envi2("NA2_PULAY_START",5);
    par.KerkerG0=envd2("NA2_KERKER",0.0);

    qchem::SolidCalcOptions o;
    o.label="Na2 AFM-seeded singlet";
    o.Nelec=2; o.multiplicity=1;                              // the explicit two-channel singlet, nUp=nDn=1
    o.species={{"Na",1}};
    o.seed=qchem::ChargeDensity::SeedStrategy::SAD;           // neutral Na: the library's 3s^1 spin pair
    o.imposeSymmetry=true;                                    // S3 resolves the SHUBNIKOV group from the flips
    o.xcMesh=qcMesh::BeckeXCParams(20, -1.0, 11);             // coarse: this gate tests the POSTCONDITION
    o.xcMesh.cellKind=qcMesh::UnitCellKind::Becke;            // PINNED -- see the header (Auto would pick Uniform)

    if (std::getenv("NA2_TRACE"))
        o.onIteration=[](const qchem::SCFIterator::SCFProgress& p)
        {
            std::cout << "[Na2 iter] " << std::setw(4) << p.iteration
                      << "  E=" << std::setprecision(10) << std::setw(15) << p.energy
                      << std::setprecision(4)
                      << "  dE=" << std::setw(11) << p.dE
                      << "  drho=" << std::setw(11) << p.drho
                      << "  [F,D]=" << std::setw(11) << p.commutator << std::endl;
        };
    // ★★★ THE MIXING STEP IS SWEPT, NOT PINNED (2026-08-27, plan step 7) -- and that is the fix for the
    // flakiness, not another lucky constant.  Deleting the pair-stream cache changed the D-aware ACTIVE SET
    // (an eps-level change, both sets eps-valid), and re-measuring alpha against it found NO PLATEAU at all:
    //
    //   alpha        0.35   0.5    0.65   0.7    0.75   0.8
    //   contracted   cap    cap    264 ✔  221 ✔  cap    248 ✔
    //   walk (=0)    cap    cap    cap    cap    cap    cap     <- ALL FOUR swept steps cap on the walk
    //
    // ⇒ WHETHER THIS FIXTURE CONVERGES IS A COIN FLIP ON THE TRUNCATION, and tuning alpha is chasing noise.
    // What is NOT a coin flip is the PHYSICS the gate is named for: in EVERY arm above -- converged or
    // capped, walk or contracted -- the moment DIES (step 23-69) and E settles on the same
    // -0.33204 state.  The instability is one slow, nearly-degenerate DENSITY mode that E cannot see: the
    // run's own fingerprint calls it "DENSITY-DEGENERATE (E settled, rho rotates -- benign)" while
    // DidConverge() calls it a failure.  ⇒ THE DURABLE FIX IS THE Δρ/N CONVERGENCE GATE
    // (doc/SCFStrategyPlan.md; item A4 on doc/OpenWork.md's anchor-moving roster), which would let this
    // state converge on its merits.  Until that lands, the honest thing is to stop betting the gate on ONE
    // draw: try the measured-good steps in order and take the first that converges.  The recipe is
    // INSTRUMENTAL (it exists so NotConverged does not pre-empt the OrderLost this gate is named for), so
    // sweeping it changes nothing about what is under test.  NA2_ALPHA pins a single value for an
    // investigation.
    // ⚠ AND BE HONEST ABOUT WHAT THE SWEEP BUYS: four draws instead of one on the DEFAULT (contracted)
    // path.  It does NOT rescue GPW_CONTRACT_CUBE=0, where all four steps hit the cap -- the reference
    // walk simply does not settle this fixture's density mode any more.  That is a property of the mode,
    // not of the kernel (the moment still dies at step 9 and E still lands on -0.332045), and it is the
    // same thing A4 would fix.
    std::vector<double> alphas{0.7, 0.8, 0.65, 0.5};
    if (const char* fixed=std::getenv("NA2_ALPHA")) alphas.assign(1, std::atof(fixed));
    std::unique_ptr<qchem::SolidCalculation> calcp;
    for (double alpha : alphas)
    {
        par.StartingRelaxRo=alpha;
        calcp=std::make_unique<qchem::SolidCalculation>(lat, MakeBasisLowQ(cell, BasisSetData::VALENCE_LOWQ_SR),
                                                       o, par);
        std::cout << "[Na2 T2] alpha="<<alpha<<" converged="<<calcp->DidConverge()
                  << " iters="<<calcp->IterationCount() << std::endl;
        if (calcp->DidConverge()) break;
    }
    const qchem::SolidCalculation& calc=*calcp;
    auto r=calc.Result();
    std::cout << "[Na2 T2] converged="<<calc.DidConverge()<<" iters="<<calc.IterationCount()
              << " E(last iterate)="<<calc.LastIterateTerms().GetTotalEnergy()
              << "\n[Na2 T2] "<<calc.Diagnostics().Summary() << std::endl;

    // THE RUN MUST CONVERGE -- otherwise this gate silently stops testing what it is named for.  With the
    // sweep above this asserts that NOT ONE of four measured-good mixing steps settled it, which is a real
    // regression signal rather than one unlucky draw.
    ASSERT_TRUE(calc.DidConverge()) << "no swept mixing step converged, so the postcondition below is "
                                       "unexercised again: " << calc.Diagnostics().Summary();
    ASSERT_FALSE(r) << "an imposed, magnetically-seeded run that relaxes to m=0 must NOT hand back an energy";
    EXPECT_EQ(r.Error().why, qchem::SCFFailure::Why::OrderLost)
        << "the failure must be the POSTCONDITION and nothing else: " << r.Error().details;

    // ...and the diagnostics must AGREE with the verdict, measured on the INTEGRATED site moment: the
    // seed carried real order (~0.47 e in the Becke basins), the answer carries none.
    const auto& diag=calc.Diagnostics();
    EXPECT_TRUE(diag.HasOrder())      << "the seed must be measurably magnetic, or the check has no teeth";
    EXPECT_GT(diag.OrderPeak(), 0.2)  << "the RAW seed's staggering, before Init's first fill eats it";
    EXPECT_TRUE(diag.OrderCollapsed());
    EXPECT_LT(diag.OrderFinal(), 0.01*diag.OrderPeak());
    // The charge channel must stay QUIET on a run that is merely non-magnetic -- the two detectors have
    // to discriminate, or the more specific one is worthless.
    EXPECT_FALSE(diag.ChargeSloshed()) << "a healthy converged run must not read as a charge runaway";
}
