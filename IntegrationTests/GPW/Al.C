// File: IntegrationTests/GPW/Al.C  Al FCC -- the degenerate-shell / Fermi-surface metal (3s2 3p1, a=7.653): smearing, annealing, global mu, IBZ.
//
// THE GRID (doc/TestSuitePlan.md §3): TEST(<Basis>_<Material>, <k>_[Grid]_[Fit]_[Sym]_[Spin]_[Occ]_[Machinery]_[Ansatz]_[Seed]_<Claim>)
// -- axis tokens in fixed order, the facade's defaults elided (k always named); the claim is CP2K (oracle anchor),
// Anchor (did-E-move), eq<Token> (a ONE-axis twin: this point equals the same point with that axis moved) or a
// property verb.  Tests are laid out in axis order.  `scripts/testgrid` renders the coverage table from these names.
//
//   GPW_Al.Γ_Imp_Stalls
//   GPW_Al.Γ_Imp_Anneal_Anchor
//   GPW_Al.k222_Smear_GlobalMu_Anchor
//   GPW_Al.k222_Imp_Smear_GlobalMu_eqFree

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


// (item 2a) THE MOTIVATION: integer aufbau CANNOT converge the degenerate 3p.  With smearing OFF, aufbau must
// place the lone 3p electron in ONE of the three degenerate p orbitals -- an arbitrary, symmetry-broken pick.
// The TOTAL ENERGY settles (|ΔE/E|~1e-13, the gap column ~0 => no frontier gap, the metallic signature) but
// the DENSITY rotates freely within the degenerate manifold, so |Δρ| floors well above tolerance and never
// converges.  This is the honest reason the smearing/annealing path below exists (mirrors the documented
// GPW_SiBox.Γ_Imp_Uni_eqFinite degenerate-shell behaviour, now for a periodic lattice).
TEST(GPW_Al, Γ_Imp_Stalls)
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
TEST(GPW_Al, Γ_Imp_Anneal_Anchor)
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
TEST(GPW_Al, k222_Smear_GlobalMu_Anchor)
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
// on the non-negative ρ_DM raster).  So the reduced run reproduces the full-mesh GPW_Al.k222_Smear_GlobalMu_Anchor free energy
// to grid/SCF tolerance (measured ~6e-8) with fewer k-points -- the IBZ payoff, done exactly (doc/GPWPlan1 item 3).
TEST(GPW_Al, k222_Imp_Smear_GlobalMu_eqFree)
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
        << "IBZ-reduced must reproduce the full-mesh free energy (GPW_Al.k222_Smear_GlobalMu_Anchor prints the same value)";
}
