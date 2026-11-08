// File: UnitTests/M_DFT.C  Molecular un-polarized Xalpha-DFT total-energy tests (facade-driven).
//
// Production path: qchem::MoleculeCalculation builds the BasisSet::Gaussian::Factory PolarizedGaussian (M&D)
// basis; the Xalpha DFT Hamiltonian adds the fitted Coulomb + Xalpha exchange on a numerical mesh.  We
// drive a full SCF for N2 and water through the facade and check the converged total energy against a
// reference -- the same total-energy check used for the atom A_DFT tests.
//
// Migrated off the QchemTester/TestMolecule scaffold (OpenWork D).  The facade auto-selects the SAD seed
// AND DIIS-from-start for DFT (CalcOptions: model is DFT -> EMax=100 so the [F,D] gate never blocks DIIS),
// which is exactly the DIIS_FromStart + SAD recipe the scaffold set by hand here -- so anchors are
// byte-for-byte the scaffold values.
//
// RE-ANCHORED TWICE on 2026-10-31 (lane B), each move judged per pin 10 against PySCF (scripts/pyscf_molecule_oracle.py):
//   1. the CONVERGED default XC mesh (MHL 80 / Lebedev 29; 30 / 5 was 20 mHa/atom low on water);
//   2. the fit bases keep their D shells (they dropped them: the Coulomb fit 4-62 mHa low) AND the Vxc default is the δ basis
//      (user ruling "δ default, A1 kept").  The INDEPENDENT ROUTE is PySCF with the SAME A1_coul Coulomb fit and exact XC: the
//      default (δ) anchors sit on it to 6-16 μHa.  The A1 Vxc fit (VxcFit::Auxiliary) keeps its own gates below, 1-6 mHa off.
#include "gtest/gtest.h"
#include <cmath>
#include <stdexcept>

import qchem.MoleculeCalculation;            // MoleculeCalculation, CalcOptions, Model, Angular
import qchem.SCFCap;                             // IterationBudget: the convergence PANEL's budget (doc/SCFRobustnessPlan.md §4)
import qchem.SCFParams;              // SCFParams (the O2 runs converge once, on their own recipe)
import qchem.StructureData;        // GetMolecule("H2O") -- the shared water geometry (D-MAKEWATER)
import qchem.Structure;              // Molecule, Atom
import qchem.Types;                  // Vector3D
import qchem.Hamiltonian.Types;   // VxcFit
using namespace qchem;

static Molecule MakeN2()
{
    Molecule m;
    m.Insert(Site{7, Vector3D<double>(-1.03, 0, 0)});
    m.Insert(Site{7, Vector3D<double>( 1.04, 0, 0)});
    return m;
}
static Molecule MakeO2()         // O2 at the experimental bond length ~1.208 Ang = 2.282 bohr, along z
{
    Molecule m;
    m.Insert(Site{8, Vector3D<double>(0, 0,  1.141)});
    m.Insert(Site{8, Vector3D<double>(0, 0, -1.141)});
    return m;
}

// Anchors are the converged Xalpha total energies (deterministic regression sentinels -- "did E
// move", NOT the physical energy; exchange-only Xalpha at these alphas is not calibrated to it).
TEST(M_DFT, N2)
{
    MoleculeCalculation calc(MoleculeSetup::Resolve(MakeN2(), {.basis = "dzvp", .model = Model::Xalpha, .xalpha = 0.75197}));
    // δ default, D shells kept: -109.228113.  Was -109.225960 (A1 Vxc fit, d-less Coulomb fit) and -109.21679 at 2e-3 RELATIVE
    // (0.2 Ha!) on the 30/5 mesh.
    EXPECT_NEAR(calc.Energy(), -109.2281129068, 1e-5);
    EXPECT_NEAR(calc.Energy(), -109.22810738, 1e-4) << "ORACLE: PySCF, the same A1_coul Coulomb fit, exact XC";
}

TEST(M_DFT, Water)
{
    MoleculeCalculation calc(MoleculeSetup::Resolve(qchem::StructureData::GetMolecule("H2O"), {.basis = "dzvp", .model = Model::Xalpha, .xalpha = 0.74000}));
    // CONVERGED Xalpha total energy (the facade's auto DIIS-from-start): [F,D] and Δρ reach ~1e-13 by
    // ~iter 14, virial V/K = -2.006 (the correct ~-2, vs the buggy -2.11).  A true converged regression
    // sentinel.  (Was -76.123348 interim limit-cycle, -79.414120 stale-cache, -76.149301 on the 30/5 mesh -- 48 mHa low, inside
    // the old 2e-3 relative -- and -76.100897 on the A1 Vxc fit with a d-less Coulomb fit.)
    EXPECT_NEAR(calc.Energy(), -76.1069880486, 1e-5);
    EXPECT_NEAR(calc.Energy(), -76.10697987, 1e-4) << "ORACLE: PySCF, the same A1_coul Coulomb fit, exact XC";
}

// --- The spherical-Gaussian (PG_Spherical) basis through the DFT path -----------------------------
// angular=Spherical selects the real-solid-harmonic orbital basis AND a spherical fit (auxiliary) basis
// (PG_Spherical::EFit_IBS, via the orbital IBS's CreateCDFit/CreateVxcFit).  The spherical Xalpha energy
// lands ~0.09 Ha ABOVE the Cartesian -- larger than HF's ~3 mHa orbital-only gap because the FIT basis
// also goes spherical (contaminant-free -> fewer auxiliary functions -> a different, poorer density fit).
// That gap is a genuine basis difference, not an error: the fit kernels are validated correct by
// M_Spherical.fit_kernels.  So this is purely a "did E move" regression anchor, with a loose sanity bound.
TEST(M_DFT, WaterSpherical)
{
    MoleculeCalculation calc(MoleculeSetup::Resolve(qchem::StructureData::GetMolecule("H2O"), {.basis   = "dzvp", .model   = Model::Xalpha,
                                   .angular = Angular::Spherical, .xalpha  = 0.74000}));
    EXPECT_NEAR(calc.Energy(), -76.1069880486, 0.05);                 // sanity: same ballpark as Cartesian
    // CONVERGED anchor (δ default, D shells kept).  PySCF spherical dzvp with the Cartesian A1_coul J-fit: -76.10758175 (16 μHa;
    // our spherical fit basis is the spherical view, so not the identical aux set).  0.6 mHa BELOW the Cartesian.
    EXPECT_NEAR(calc.Energy(), -76.1075980849, 1e-5);
}

// Polarized Xalpha water through the facade's DEFAULT (auto SAD) seed.  This used to SEGFAULT: the SAD seed
// is a spin-agnostic total density, and the polarized Vxc null-derefed when its dynamic_cast<Polarized_CD>
// of the seed failed (the assert was compiled out in Release).  Fixed in FittedVxcPol::GetMatrix -- a
// spin-unpolarized seed (rho_up=rho_down=rho/2) maps each channel's iteration-0 Vxc to the unpolarized Vxc
// of the total density.  Water is closed shell, so with the SAME alpha the polarized run must converge to
// the SAME energy as the unpolarized anchor -- and it does, to ~1e-11 (confirming a correct seed Fock).
TEST(M_DFT, WaterPolarizedSAD)
{
    MoleculeCalculation calc(MoleculeSetup::Resolve(qchem::StructureData::GetMolecule("H2O"), {.basis = "dzvp", .model = Model::Xalpha, .spin = SpinGroup::Polarized, .xalpha = 0.74000}));
    // Polarized water Xalpha is the oscillatory case (commutator ~2-4; see M_Sym RunDFT) -- it needs more
    // iterations + damping than the facade's quick default, so re-Converge with tight params (DIIS already
    // drives from the start for DFT).
    ASSERT_TRUE(calc.Converge({.NMaxIter = 60, .MinΔρPerElectron = 1e-7, .MinΔFD = 1e-9, .MinVirial = 1e2, .MinFD = 1e-7,
                   .StartingRelaxRo = 0.5, .MergeTol = 1e-4, .Verbose = false}));
    EXPECT_NEAR(calc.Energy(), -76.1069880486, 1e-5);   // == the unpolarized Xalpha anchor
}

// Spin-native (polarized) parameter-free LDA water (Dirac exchange + spin-native VWN5 correlation) through
// the facade's auto SAD seed -- OpenWork B2: Ham_DFTcorr_P + FittedVcorrPol.  Water is closed shell, so at
// convergence zeta=0 everywhere and the polarized run must collapse to the UNPOLARIZED LDA anchor
// (M_Calculation.WaterLDA = -75.8794551).  This exercises the new polarized correlation term end to end:
// the SAD-seed fallback (rho_up=rho_down=rho/2), the per-spin v_c^sigma fit, and the two-channel E_c.
TEST(M_DFT, WaterPolarizedLDA)
{
    MoleculeCalculation calc(MoleculeSetup::Resolve(qchem::StructureData::GetMolecule("H2O"), {.basis = "dzvp", .model = Model::LDA, .spin = SpinGroup::Polarized}));
    ASSERT_TRUE(calc.Converge({.NMaxIter = 60, .MinΔρPerElectron = 1e-7, .MinΔFD = 1e-9, .MinVirial = 1e2, .MinFD = 1e-7,
                   .StartingRelaxRo = 0.5, .MergeTol = 1e-4, .Verbose = false}));
    EXPECT_NEAR(calc.Energy(), -75.8794551388, 1e-6);   // collapses to the unpolarized LDA anchor (closed shell)
}

// OPEN-SHELL through the facade multiplicity knob (OpenWork B4): triplet O2 (16 e-, multiplicity 3 =>
// nUp=9, nDown=7).  multiplicity>1 auto-promotes the calc to SpinGroup::Polarized, builds Molecule_EC(9,7), and
// runs spin-native LDA.  The spin texture is real here (the half-filled pi* gives the triplet), so this
// genuinely exercises zeta!=0 -- not the closed-shell collapse the WaterPolarized tests hit.  "Did E move"
// regression sentinel (converged spin-native LDA total energy).
TEST(M_DFT, OxygenTripletLDA)
{
    MoleculeCalculation calc(MoleculeSetup::Resolve(MakeO2(), {.basis = "dzvp", .model = Model::LDA, .multiplicity = 3}));
    ASSERT_TRUE(calc.Converge({.NMaxIter = 80, .MinΔρPerElectron = 1e-7, .MinΔFD = 1e-9, .MinVirial = 1e2, .MinFD = 1e-7,
                   .StartingRelaxRo = 0.5, .MergeTol = 1e-4, .Verbose = false}));
    // Converged spin-native LDA triplet (δ default, D shells kept).  Was -149.263716 on the A1 per-spin Vxc fit with a d-less
    // Coulomb fit -- 35.5 mHa above PySCF exact, which a row once blamed on the polarized machinery; it was the fit.
    EXPECT_NEAR(calc.Energy(), -149.3006843396, 1e-4);
    EXPECT_NEAR(calc.Energy(), -149.30067460, 1e-4) << "ORACLE: PySCF UKS, the same A1_coul Coulomb fit, exact XC";
}

// REGRESSION for the O2-HF-triplet SEGV: the polarized level display (tCompositeWF::DisplayEigenPolarized) crashed on
// an open-shell level present in only ONE spin channel (find() is UB-in-Release on a miss).  Fixed via
// EnergyLevels::FindOrNull + labelling from the always-valid combined el.qns.  This runs it with Verbose ON
// (the crash path) and pins the converged HF triplet energy.  "Did E move" + no-crash.
TEST(M_DFT, OxygenTripletHF_VerboseDisplay)
{
    MoleculeCalculation calc(MoleculeSetup::Resolve(MakeO2(), {.basis = "dzvp", .model = Model::HF, .multiplicity = 3}));
    ASSERT_TRUE(calc.Converge({.NMaxIter = 80, .MinΔρPerElectron = 1e-6, .MinΔFD = 1e-8, .MinVirial = 1e2, .MinFD = 1e-6,
                   .StartingRelaxRo = 0.5, .MergeTol = 1e-4, .Verbose = true}));   // Verbose => DisplayEigen (the crash path)
    EXPECT_NEAR(calc.Energy(), -149.6295470276, 1e-4);   // converged HF triplet O2
}

// Physics check that the spin channel is doing real work: the triplet ground state must lie BELOW the
// closed-shell singlet (Hund's rule).  A flat/zero gap would mean the open-shell occupation silently
// collapsed to closed shell -- the failure mode B3/B4 exist to prevent.
TEST(M_DFT, OxygenTripletBelowSinglet_Panel)   // PANEL case 11: a degenerate partly filled level (smeared singlet)
{
    qchem::IterationBudget budget(55, "O2 triplet + smeared singlet (molecular, degenerate pi*; 46 measured at kT 1e-2, 2026-10-31 -- "
                                      "was 26 at kT 1e-4 with the d-less Coulomb fit)");
    // The closed-shell SINGLET has two electrons for the exactly degenerate pi*x / pi*y pair: a single integer determinant
    // flips them between the two every iteration and never converged (2026-10-29, the cap guard).  A smearing occupies the
    // pair EVENLY -- the symmetric closed-shell singlet, still above the triplet by the exchange the triplet gains.  kT 1e-2,
    // NOT a fudge (2026-10-31): an exactly degenerate open shell has no integer-aufbau fixed point, and kT must exceed the
    // splitting the anisotropic Hartree re-opens -- kT 1e-4 converged only while the fit basis lacked its D shells (a spherical
    // Hartree); with them it sloshes at kT 1e-4 and 1e-3 under DIIS, GDM and Ladder alike ("smear wider than the splitting you
    // are curing", GPW_SiBox).  Compared on the INTERNAL energy E = A - (-TS), the free energy's entropy left out.
    // Integer occupations for the triplet (pi*x up, pi*y up: no tie).
    auto run = [](int mult) {
        SCFParams p{.NMaxIter = 100, .MinΔρPerElectron = 1e-7, .MinVirial = 1e2, .MinFD = 1e-7,
                    .StartingRelaxRo = 0.5, .MergeTol = 1e-4, .Verbose = false};
        if (mult==1) p.SmearingkT=1e-2;   // open row: doc/OpenWork.md §3 "Molecular SCF cannot converge a DEGENERATE open-shell density at small kT"
        MoleculeCalculation c(MoleculeSetup::Resolve(MakeO2(), {.basis = "dzvp", .model = Model::LDA, .multiplicity = mult}), p);
        EXPECT_TRUE(c.IsConverged()) << "multiplicity " << mult;
        const auto E=c.EnergyTerms();
        return E.GetTotalEnergy()-E["MinusTS"];   // the internal energy (the -TS of a smeared run left out; 0 for the triplet)
    };
    const double triplet = run(3);   // nUp=9, nDown=7  (SpinGroup::Polarized, auto)
    const double singlet = run(1);   // nUp=8, nDown=8  (closed shell)
    EXPECT_LT(triplet, singlet - 0.01);   // triplet is the ground state, by a clear margin
}

// A multiplicity whose parity disagrees with Ne must be rejected, not silently miscounted: water has 10
// electrons (even), so a doublet (multiplicity 2 => 2S=1, odd) is impossible.  Fail loud (B4 validation).
TEST(M_DFT, BadMultiplicityThrows)
{
    EXPECT_THROW(MoleculeCalculation(MoleculeSetup::Resolve(qchem::StructureData::GetMolecule("H2O"), {.basis = "sto-3g", .model = Model::LDA, .multiplicity = 2})),
                 std::runtime_error);
}

// THE AUXILIARY Vxc FIT (VxcFit::Auxiliary -- the Gaussian A1 set), kept as a strategy since δ became the default (user ruling
// "δ default, A1 kept"; pin 30: it keeps the ONE XC term honest over every fit basis).  The same FittedVxc consumes it through
// the fit face.  Its error is the A1 fit of v_xc AND eps_xc: 1-6 mHa against PySCF with the SAME A1_coul Coulomb fit and exact XC.
static Molecule MakeF2()   // ~/Code/cp2k-runs/F2 geometry (bohr)
{
    Molecule m;
    m.Insert(Site{9, Vector3D<double>(0, 0, -1.3343)});
    m.Insert(Site{9, Vector3D<double>(0, 0,  1.3343)});
    return m;
}

TEST(M_DFT, WaterLDA_AuxiliaryFit)
{
    MoleculeCalculation calc(MoleculeSetup::Resolve(qchem::StructureData::GetMolecule("H2O"), {.basis = "dzvp", .model = Model::LDA,
                                                                                       .vxcFit = Hamiltonian::VxcFit::Auxiliary}));
    ASSERT_TRUE(calc.Result());
    EXPECT_NEAR(calc.Energy(), -75.8808368818, 1e-5);   // anchor
    EXPECT_NEAR(calc.Energy(), -75.87944736, 3e-3) << "ORACLE: PySCF A1_coul J-fit route; the gap (-1.4 mHa) is the A1 Vxc fit";
}

TEST(M_DFT, OxygenTripletLDA_AuxiliaryFit)
{
    MoleculeCalculation calc(MoleculeSetup::Resolve(MakeO2(), {.basis = "dzvp", .model = Model::LDA, .multiplicity = 3,
                                                       .vxcFit = Hamiltonian::VxcFit::Auxiliary}));
    ASSERT_TRUE(calc.Converge({.NMaxIter = 80, .MinΔρPerElectron = 1e-7, .MinΔFD = 1e-9, .MinVirial = 1e2, .MinFD = 1e-7,
                               .StartingRelaxRo = 0.5, .MergeTol = 1e-4, .Verbose = false}));
    EXPECT_NEAR(calc.Energy(), -149.3018824495, 1e-4);   // anchor
    EXPECT_NEAR(calc.Energy(), -149.30067460, 3e-3) << "ORACLE: PySCF UKS A1_coul J-fit route; the gap (-1.2 mHa) is the A1 Vxc fit";
}

// F2 GTH q7 (the case the "Molecular PP route is 4 mHa/atom off CP2K" row was filed on) on BOTH routes: the δ default is CP2K's
// number to 1.4 mHa (PySCF's own A1_coul J-fit route to 57 μHa) -- the PP route was never the defect.
TEST(M_DFT, F2_PP)
{
    MoleculeCalculation calc(MoleculeSetup::Resolve(MakeF2(), {.basis = "valence_lowq_vb", .model = Model::LDA, .pp = "gth-lda-valence"}));
    ASSERT_TRUE(calc.Result());
    // Re-judged 2026-11-05 (pin 10): the closed-form PP integrals moved this +4.2e-5, from 5.7e-5 to 1.5e-5 off PySCF.
    EXPECT_NEAR(calc.Energy(), -48.2018945557, 1e-5);   // anchor (δ default)
    EXPECT_NEAR(calc.Energy(), -48.20187971, 1e-4) << "ORACLE: PySCF, the same A1_coul Coulomb fit, exact XC (CP2K -48.20058)";
}
TEST(M_DFT, F2_PP_AuxiliaryFit)
{
    MoleculeCalculation calc(MoleculeSetup::Resolve(MakeF2(), {.basis = "valence_lowq_vb", .model = Model::LDA, .vxcFit = Hamiltonian::VxcFit::Auxiliary,
                                                       .pp = "gth-lda-valence"}));
    ASSERT_TRUE(calc.Result());
    EXPECT_NEAR(calc.Energy(), -48.1970701902, 1e-5);   // anchor (re-judged with F2_PP above, the same +4.1e-5)
    EXPECT_NEAR(calc.Energy(), -48.20187971, 7e-3) << "ORACLE: PySCF A1_coul J-fit route; the gap (+4.8 mHa) is the A1 Vxc fit";
}

// THE δ HARTREE ROUTE (CoulombFit::Delta): ρ on the XC mesh, the direct Coulomb sum, first-order energy, through the same FittedVee.
// The direct sum is coarse, by measurement (TeF₆ prototype: J to about 1e-2), so the oracle tolerance says so.
TEST(M_DFT, WaterLDA_DeltaCoulomb_Long)
{
    MoleculeCalculation calc(MoleculeSetup::Resolve(qchem::StructureData::GetMolecule("H2O"), {.basis = "dzvp", .model = Model::LDA,
                                                                                       .coulomb = {.fit = Hamiltonian::CoulombFit::Delta, .solver = Hamiltonian::CoulombSolver::DirectSum}}));
    ASSERT_TRUE(calc.Result());
    // The mesh pilot sizes XC and leaves this Hartree out (its error does not fall with the mesh, HartreeConvergesWithMesh): chasing it
    // walked the ladder at O(N²) per rung and ran ~10 h (2026-11-08).  A settled climb is the proof it was not chased.
    EXPECT_TRUE(calc.PilotCost().settled);
    EXPECT_NEAR(calc.Energy(), -75.8917526160, 1e-5);   // anchor
    EXPECT_NEAR(calc.Energy(), -75.87729931, 2e-2) << "ORACLE: PySCF, exact J and XC; the gap (-14.5 mHa) is the direct sum";
}

// ...and with the radial Poisson solve, the route to use: the same mesh, the operator at the mesh's own limit.  It sits 15x closer
// to exact J than the Gaussian A1 fit does (-2.1 mHa, the PySCF A1_coul row).  RE-JUDGED 2026-11-08 (pin 10, lane B): the unstated mesh
// now climbs to eps (GeometryOptimizationPlan §14): δ-ρ Hartree rides it and misses eps x N at (80, 29), so the pilot raises it to
// (100, 29), and the energy moves -39.9 μHa, from +30.6 to -9.3 μHa off the exact-J oracle.
TEST(M_DFT, WaterLDA_DeltaCoulomb_Poisson)
{
    MoleculeCalculation calc(MoleculeSetup::Resolve(qchem::StructureData::GetMolecule("H2O"), {.basis = "dzvp", .model = Model::LDA,
                                                                                       .coulomb = {.fit = Hamiltonian::CoulombFit::Delta}}));
    ASSERT_TRUE(calc.Result());
    EXPECT_NEAR(calc.Energy(), -75.8773086479, 1e-6);   // anchor (re-judged 2026-11-08, above)
    EXPECT_NEAR(calc.Energy(), -75.87729931, 5e-5) << "ORACLE: PySCF, exact J and XC; the gap (-9 μHa) is the mesh's, inside eps x N";
}

// The multipole fit has δ x Poisson's potential, so its SCF is the same one; only the energy's second leg differs.  Re-judged with the
// δ case above (2026-11-08): the pilot's (100, 29) moves it -39.5 μHa, from +32 to -7.9 μHa off the oracle.
TEST(M_DFT, WaterLDA_MultipoleCoulomb)
{
    MoleculeCalculation calc(MoleculeSetup::Resolve(qchem::StructureData::GetMolecule("H2O"), {.basis = "dzvp", .model = Model::LDA,
                                                                                       .coulomb = {.fit = Hamiltonian::CoulombFit::Multipole}}));
    ASSERT_TRUE(calc.Result());
    EXPECT_NEAR(calc.Energy(), -75.8773072275, 1e-6);   // anchor (re-judged 2026-11-08, above)
    EXPECT_NEAR(calc.Energy(), -75.87729931, 5e-5) << "ORACLE: PySCF, exact J and XC: delta's level, -8 μHa.  This rung was 530 out when rho-tilde was synthesized onto the mesh";
}

// The multipole route's energy depends less on the l cut than δ's: at 150/35, 11 μHa between l = 12 and 16 where δ moves 35.
TEST(M_DFT, WaterLDA_MultipoleCoulomb_eqLCut_Long)
{
    auto energy=[](int lMax)
    {
        MoleculeCalculation calc(MoleculeSetup::Resolve(qchem::StructureData::GetMolecule("H2O"),
            {.basis = "dzvp", .model = Model::LDA, .mesh = {.radial = qcMesh::RadialKind::MHL, .nRadial = 150, .mhl_m = 3, .mhl_alpha = 2.0,
                                                      .angular = qcMesh::AngularKind::Lebedev, .angularDegree = 35, .beckeOrder = 2},
             .coulomb = {.fit = Hamiltonian::CoulombFit::Multipole, .lMax = lMax}}));
        EXPECT_TRUE(calc.Result());
        return calc.Energy();
    };
    const double e12=energy(12);
    EXPECT_NEAR(e12, energy(16), 2e-5);
    EXPECT_NEAR(e12, -75.87729931, 2e-5) << "ORACLE: PySCF, exact J and XC";
}

// Multipole for Vxc, through the fitter as the A1 route is.  The Hartree stays on the default route, so only the XC term changed.
TEST(M_DFT, WaterLDA_MultipoleVxc)
{
    MoleculeCalculation calc(MoleculeSetup::Resolve(qchem::StructureData::GetMolecule("H2O"), {.basis = "dzvp", .model = Model::LDA,
                                                                                       .vxcFit = Hamiltonian::VxcFit::Multipole}));
    ASSERT_TRUE(calc.Result());
    EXPECT_NEAR(calc.Energy(), -75.87944736, 1e-2) << "the A1 Hartree value; the gap is the multipole Vxc's";
}
