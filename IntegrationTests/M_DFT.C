// File: UnitTests/M_DFT.C  Molecular un-polarized Xalpha-DFT total-energy tests (facade-driven).
//
// Production path: qchem::MoleculeCalculation builds the BasisSet::Gaussian::Factory PolarizedGaussian (M&D)
// basis; the Xalpha DFT Hamiltonian adds the fitted Coulomb + Xalpha exchange on a numerical mesh.  We
// drive a full SCF for N2 and water through the facade and check the converged total energy against a
// reference -- the same total-energy check used for the atom A_DFT tests.
//
// Migrated off the QchemTester/TestMolecule scaffold (OpenWork D).  The facade auto-selects the SAD seed
// AND DIIS-from-start for DFT (RunOptions: model is DFT -> EMax=100 so the [F,D] gate never blocks DIIS),
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
#include <iostream>
#include <iomanip>
#include <string>

import qchem.MoleculeCalculation;            // MoleculeCalculation, RunOptions, Model, Angular
import qchem.SCFCap;                             // IterationBudget: the convergence PANEL's budget (doc/SCFRobustnessPlan.md §4)
import qchem.SCFParams;              // SCFParams (the O2 runs converge once, on their own recipe)
import qchem.StructureData;        // GetMolecule("H2O") -- the shared water geometry (D-MAKEWATER)
import qchem.Structure;              // Molecule, Atom
import qchem.Types;                  // Vector3D
import qchem.Hamiltonian.Types;   // VxcFit
import qchem.ChargeDensity.Seed;   // SeedStrategy (the Ionic SAD gates)
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
    auto calc=MoleculeCalculation::Prepare(MoleculeSetup::Resolve(MakeN2(), {.model=Model::Xalpha, .xalpha=0.75197, .molecule={.basis="dzvp"}})).TakeValue(); (void)calc->Converge();
    // RE-JUDGED 2026-11-22 (pin 10, lane S; user ruling "Yes, δ default"): the Hartree is the δ fit now, not A1_coul.  From -109.2281129068
    // (the A1_coul Coulomb fit, 0.92 mHa low); PySCF EXACT J and XC -109.22719891, 6.1 μHa away.  Earlier: -109.225960 (A1 Vxc fit, d-less
    // Coulomb fit) and -109.21679 at 2e-3 RELATIVE (0.2 Ha!) on the 30/5 mesh.
    EXPECT_NEAR(calc->Energy(), -109.2271928297, 1e-5);
    EXPECT_NEAR(calc->Energy(), -109.22719891, 2e-5) << "ORACLE: PySCF, exact J and XC";
}

TEST(M_DFT, Water)
{
    auto calc=MoleculeCalculation::Prepare(MoleculeSetup::Resolve(qchem::StructureData::GetMolecule("H2O"), {.model=Model::Xalpha, .xalpha=0.74000, .molecule={.basis="dzvp"}})).TakeValue(); (void)calc->Converge();
    // CONVERGED Xalpha total energy (the facade's auto DIIS-from-start): [F,D] and Δρ reach ~1e-13 by
    // ~iter 14, virial V/K = -2.006 (the correct ~-2, vs the buggy -2.11).  A true converged regression
    // sentinel.  (Was -76.123348 interim limit-cycle, -79.414120 stale-cache, -76.149301 on the 30/5 mesh -- 48 mHa low, inside
    // the old 2e-3 relative -- and -76.100897 on the A1 Vxc fit with a d-less Coulomb fit.)
    // RE-JUDGED 2026-11-22 (pin 10, lane S; user ruling "Yes, δ default"): from -76.1069880486 (the A1_coul Coulomb fit, 2.2 mHa low);
    // PySCF EXACT J and XC -76.10476847, 9.6 μHa away.
    EXPECT_NEAR(calc->Energy(), -76.1047780883, 1e-5);
    EXPECT_NEAR(calc->Energy(), -76.10476847, 3e-5) << "ORACLE: PySCF, exact J and XC";
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
    auto calc=MoleculeCalculation::Prepare(MoleculeSetup::Resolve(qchem::StructureData::GetMolecule("H2O"), {.model=Model::Xalpha, .xalpha=0.74000, .molecule={.basis="dzvp", .angular=Angular::Spherical}})).TakeValue(); (void)calc->Converge();
    EXPECT_NEAR(calc->Energy(), -76.1069880486, 0.05);                 // sanity: same ballpark as Cartesian
    // RE-JUDGED 2026-11-22 (pin 10, lane S; user ruling "Yes, δ default"): the Hartree is the δ fit, so no spherical view of a Cartesian fit
    // file is involved any more.  From -76.1075980849 (the spherical A1_coul view, 3.5 mHa low); PySCF spherical dzvp, EXACT J and XC
    // -76.10410443, 9.8 μHa away.  0.66 mHa ABOVE the Cartesian (the spherical orbital space is the smaller one).
    EXPECT_NEAR(calc->Energy(), -76.1041141992, 1e-5);
    EXPECT_NEAR(calc->Energy(), -76.10410443, 3e-5) << "ORACLE: PySCF spherical, exact J and XC";
}

// THE TIGHT RE-CONVERGE these tests run after the constructor's default one.  It lives in the SETUP (\c RunOptions::scf, R10) and
// \c Converge() runs it: a whole SCFParams handed to Converge would reset every unnamed field to its library default.
static const SCFParams kTight60 = {.NMaxIter = 60, .MinΔρPerElectron = 1e-7, .MinΔFD = 1e-9, .MinVirial = 1e2, .MinFD = 1e-7, .StartingRelaxRo = 0.5, .MergeTol = 1e-4, .Verbose = false};
static const SCFParams kTight80 = {.NMaxIter = 80, .MinΔρPerElectron = 1e-7, .MinΔFD = 1e-9, .MinVirial = 1e2, .MinFD = 1e-7, .StartingRelaxRo = 0.5, .MergeTol = 1e-4, .Verbose = false};
static const SCFParams kTightHF_Verbose = {.NMaxIter = 80, .MinΔρPerElectron = 1e-6, .MinΔFD = 1e-8, .MinVirial = 1e2, .MinFD = 1e-6, .StartingRelaxRo = 0.5, .MergeTol = 1e-4, .Verbose = true};   // Verbose => DisplayEigen

// Polarized Xalpha water through the facade's DEFAULT (auto SAD) seed.  This used to SEGFAULT: the SAD seed
// is a spin-agnostic total density, and the polarized Vxc null-derefed when its dynamic_cast<Polarized_CD>
// of the seed failed (the assert was compiled out in Release).  Fixed in FittedVxcPol::GetMatrix -- a
// spin-unpolarized seed (rho_up=rho_down=rho/2) maps each channel's iteration-0 Vxc to the unpolarized Vxc
// of the total density.  Water is closed shell, so with the SAME alpha the polarized run must converge to
// the SAME energy as the unpolarized anchor -- and it does, to ~1e-11 (confirming a correct seed Fock).
TEST(M_DFT, WaterPolarizedSAD)
{
    auto calc=MoleculeCalculation::Prepare(MoleculeSetup::Resolve(qchem::StructureData::GetMolecule("H2O"), {.spin=SpinGroup::Polarized, .model=Model::Xalpha, .xalpha=0.74000, .scf=kTight60, .molecule={.basis="dzvp"}})).TakeValue(); (void)calc->Converge();
    // Polarized water Xalpha is the oscillatory case (commutator ~2-4; see M_Sym RunDFT) -- it needs more
    // iterations + damping than the facade's quick default, so re-Converge with tight params (DIIS already
    // drives from the start for DFT).
    ASSERT_TRUE(calc->Converge());
    EXPECT_NEAR(calc->Energy(), -76.1047780883, 1e-5);   // == the unpolarized Xalpha anchor (re-judged 2026-11-22, δ default)
}

// Spin-native (polarized) parameter-free LDA water (Dirac exchange + spin-native VWN5 correlation) through
// the facade's auto SAD seed -- OpenWork B2: Ham_DFTcorr_P + FittedVcorrPol.  Water is closed shell, so at
// convergence zeta=0 everywhere and the polarized run must collapse to the UNPOLARIZED LDA anchor
// (M_Calculation.WaterLDA = -75.8794551).  This exercises the new polarized correlation term end to end:
// the SAD-seed fallback (rho_up=rho_down=rho/2), the per-spin v_c^sigma fit, and the two-channel E_c.
TEST(M_DFT, WaterPolarizedLDA)
{
    auto calc=MoleculeCalculation::Prepare(MoleculeSetup::Resolve(qchem::StructureData::GetMolecule("H2O"), {.spin=SpinGroup::Polarized, .model=Model::LDA, .scf=kTight60, .molecule={.basis="dzvp"}})).TakeValue();
    ASSERT_TRUE(calc->Converge());
    EXPECT_NEAR(calc->Energy(), -75.8773086479, 1e-6);   // collapses to the unpolarized LDA anchor (closed shell; re-judged 2026-11-22, δ default)
}

// OPEN-SHELL through the facade multiplicity knob (OpenWork B4): triplet O2 (16 e-, multiplicity 3 =>
// nUp=9, nDown=7).  multiplicity>1 auto-promotes the calc to SpinGroup::Polarized, builds Molecule_EC(9,7), and
// runs spin-native LDA.  The spin texture is real here (the half-filled pi* gives the triplet), so this
// genuinely exercises zeta!=0 -- not the closed-shell collapse the WaterPolarized tests hit.  "Did E move"
// regression sentinel (converged spin-native LDA total energy).
TEST(M_DFT, OxygenTripletLDA)
{
    auto calc=MoleculeCalculation::Prepare(MoleculeSetup::Resolve(MakeO2(), {.multiplicity=3, .model=Model::LDA, .scf=kTight80, .molecule={.basis="dzvp"}})).TakeValue();
    ASSERT_TRUE(calc->Converge());
    // RE-JUDGED 2026-11-22 (pin 10, lane S; user ruling "Yes, δ default"): from -149.3006843396 (the A1_coul Coulomb fit, 1.5 mHa low);
    // PySCF UKS EXACT J and XC -149.29917753, 20 μHa away (eps x N_atoms = 2e-5).  Earlier -149.263716 on the A1 per-spin Vxc fit with a
    // d-less Coulomb fit -- 35.5 mHa above PySCF exact, which a row once blamed on the polarized machinery; it was the fit.
    EXPECT_NEAR(calc->Energy(), -149.2991575707, 1e-4);
    EXPECT_NEAR(calc->Energy(), -149.29917753, 3e-5) << "ORACLE: PySCF UKS, exact J and XC";
}

// THE SPIN-RESOLVED MOLECULAR SAD SEED on the O2 triplet (OpenWork §2 "Molecular spin-resolved SAD seed", lane B 2026-11-16).  The O entry of
// each seed library carries its Hund pair (moment 2): the valence O q6 under the GTH pseudopotential, the all-electron O since the light
// open-shell neutrals were regenerated with their pairs (the same day).  So the polarized run's Auto seed (SAD) plants each O's triplet, scaled
// by one factor to the run's moment nUp - nDown = 2 (two O triplets plant 4).  The claim: it reaches the triplet the core guess reaches.
static void SpinResolvedSeedReachesTheCoreGuessState(const RunOptions& base, const char* what)
{
    double E[2]={0,0};
    const ChargeDensity::SeedStrategy seeds[2]={ChargeDensity::SeedStrategy::Default, ChargeDensity::SeedStrategy::CoreGuess};
    for (int i=0;i<2;i++)
    {
        RunOptions o=base;
        o.seed=seeds[i];
        o.scf=kTight80;
        auto calc=MoleculeCalculation::Prepare(MoleculeSetup::Resolve(MakeO2(), o)).TakeValue();
        ASSERT_TRUE(calc->Converge()) << what << (i ? ": core guess" : ": spin-resolved SAD");
        E[i]=calc->Energy();
        std::cout << "[O2 triplet, " << what << "] " << (i ? "core guess" : "spin-resolved SAD") << ": E " << std::setprecision(12) << E[i]
                  << ", " << calc->IterationCount() << " iterations" << std::endl;
    }
    EXPECT_NEAR(E[0], E[1], 1e-7) << what << ": the spin-resolved seed reaches the core guess's triplet";
}
TEST(M_DFT, OxygenTripletPP_TheSpinResolvedSeedReachesTheCoreGuessState)
{
    SpinResolvedSeedReachesTheCoreGuessState({.multiplicity=3, .model=Model::LDA, .pp="gth-lda-valence", .molecule={.basis="valence_lowq_vb"}}, "PP");
}
TEST(M_DFT, OxygenTripletLDA_TheSpinResolvedSeedReachesTheCoreGuessState)
{
    SpinResolvedSeedReachesTheCoreGuessState({.multiplicity=3, .model=Model::LDA, .molecule={.basis="dzvp"}}, "all-electron");
}

// REGRESSION for the O2-HF-triplet SEGV: the polarized level display (tCompositeWF::DisplayEigenPolarized) crashed on
// an open-shell level present in only ONE spin channel (find() is UB-in-Release on a miss).  Fixed via
// EnergyLevels::FindOrNull + labelling from the always-valid combined el.qns.  This runs it with Verbose ON
// (the crash path) and pins the converged HF triplet energy.  "Did E move" + no-crash.
TEST(M_DFT, OxygenTripletHF_VerboseDisplay)
{
    auto calc=MoleculeCalculation::Prepare(MoleculeSetup::Resolve(MakeO2(), {.multiplicity=3, .model=Model::HF, .scf=kTightHF_Verbose, .molecule={.basis="dzvp"}})).TakeValue();
    ASSERT_TRUE(calc->Converge());   // Verbose => DisplayEigen (the crash path)
    EXPECT_NEAR(calc->Energy(), -149.6295470276, 1e-4);   // converged HF triplet O2
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
        auto c=MoleculeCalculation::Prepare(MoleculeSetup::Resolve(MakeO2(), {.multiplicity=mult, .model=Model::LDA, .scf=p, .molecule={.basis="dzvp"}})).TakeValue(); (void)c->Converge();
        EXPECT_TRUE(c->IsConverged()) << "multiplicity " << mult;
        const auto E=c->EnergyTerms();
        return E.GetTotalEnergy()-E["MinusTS"];   // the internal energy (the -TS of a smeared run left out; 0 for the triplet)
    };
    const double triplet = run(3);   // nUp=9, nDown=7  (SpinGroup::Polarized, auto)
    const double singlet = run(1);   // nUp=8, nDown=8  (closed shell)
    EXPECT_LT(triplet, singlet - 0.01);   // triplet is the ground state, by a clear margin
}

// A multiplicity whose parity disagrees with Ne must be rejected, not silently miscounted: water has 10
// electrons (even), so a doublet (multiplicity 2 => 2S=1, odd) is impossible.  Refused at the setup's user door (B4 validation).
// (Until 2026-11-14 this asked for basis "sto-3g", which is not registered, so it passed on the BASIS lookup's throw and never reached
// the multiplicity check.)
TEST(M_DFT, BadMultiplicityIsRefused)
{
    const auto refused=MoleculeSetup::Admit(qchem::StructureData::GetMolecule("H2O"), {.multiplicity=2, .model=Model::LDA, .molecule={.basis="dzvp"}});
    ASSERT_FALSE(refused);
    EXPECT_EQ(refused.Error().policy, "SpinMultiplicity");
}

// THE AUXILIARY Vxc FIT (VxcFit::Auxiliary -- the Gaussian A1 set), kept as a strategy since δ became the default (user ruling
// "δ default, A1 kept"; pin 30: it keeps the ONE XC term honest over every fit basis).  At xc.rho exact (the W cell) the one XC term
// runs it as the fitted-potential sampler: E is the exact XC of rho[D], H the A1 fit of v.  Until FittedVxc's retirement (2026-11-15)
// eps was fitted for E as well, 1-6 mHa off PySCF with the SAME A1_coul Coulomb fit and exact XC; now 2-80 μHa, second order.
static Molecule MakeF2()   // ~/Code/cp2k-runs/F2 geometry (bohr)
{
    Molecule m;
    m.Insert(Site{9, Vector3D<double>(0, 0, -1.3343)});
    m.Insert(Site{9, Vector3D<double>(0, 0,  1.3343)});
    return m;
}

// TeF6 SPHERICAL under a pseudopotential -- the case the OpenWork row "a spherical molecular PP run dies on its first H at radial
// meshes of 50 points or more" was filed on (lane B, 2026-11-05: "Invalid setup of symmetric matrix" at the first diagonalisation for
// (50,17) (50,29) (80,29) (150,29) (200,35); (30..40,17) ran, Cartesian ran every rung).  FIXED AT HEAD, CAUSE UNIDENTIFIED (lane B,
// 2026-11-12, doc/Records/OpenWork_History5.md): every rung converges in 11 iterations, and no spherical orbital or A1 fit function is
// non-finite at any mesh point or at a nucleus.  This is the crash recipe itself: the core seed, DIIS, the first failing mesh, stated.
static Molecule MakeTeF6()   // octahedron, Te-F 1.82 Angstrom (the 2026-11-05 runs' geometry)
{
    const double d=1.82/0.529177210903;
    Molecule m;
    m.Insert(Site{52, Vector3D<double>(0,0,0)});
    for (const Vector3D<double>& r : {Vector3D<double>(d,0,0), Vector3D<double>(-d,0,0), Vector3D<double>(0,d,0), Vector3D<double>(0,-d,0),
                                      Vector3D<double>(0,0,d), Vector3D<double>(0,0,-d)})
        m.Insert(Site{9, r});
    return m;
}
TEST(M_DFT, TeF6_PP_Spherical_Converges_Long)
{
    RunOptions o{.model=Model::LDA, .pp="gth-lda-valence", .molecule={.basis="dzvp_molopt_sr_gth"}};
    o.species = {{"Te", 6}};                                    // the basis's Te block is q6 (its PP-PROVENANCE line)
    o.molecule.angular    = Angular::Spherical;
    o.seed       = ChargeDensity::SeedStrategy::CoreGuess;
    MoleculeXCMesh(o).nRadial = 50; MoleculeXCMesh(o).angularDegree = 17;                // the first rung that died
    o.coulomb.auxBasis = "A1_coul";   // the recipe's Hartree, STATED: the δ default (user ruling 2026-11-22) on this coarse mesh is 35 mHa off
    auto calc=MoleculeCalculation::Prepare(MoleculeSetup::Resolve(MakeTeF6(), o)).TakeValue(); (void)calc->Converge();
    ASSERT_TRUE(calc->Result()) << "the 2026-11-05 crash recipe must converge";
    EXPECT_NEAR(calc->Energy(), -153.4149661057, 1e-5);             // anchor at the stated (50,17)
    EXPECT_NEAR(calc->Energy(), -153.4140873524, 2e-3) << "the (200,35) value: (50,17) is 0.9 mHa off it, the coarse mesh's own error";
}

// AN ALL-ELECTRON MOLECULE PAST NEON STARTS FROM SAD (OpenWork §3 "The default SAD seed stops at Ne", lane B 2026-11-13): the all-electron
// seed library held Z 1-10 only, so these took the core guess with a WARN (before 3afed2228 they died on an uncaught throw).  The library
// now holds every element NIST tabulates (scripts/gen-ae-seeds), and the Auto seed reads it: the SAD start reaches the SAME converged
// state as the core guess -- a seed changes the path, never the answer.
static Molecule MakeH2S()   // C2v, S-H 1.336 Angstrom, H-S-H 92.1 degrees (experimental)
{
    const double r=1.336/0.529177210903, h=0.5*92.1*3.14159265358979323846/180.0;
    Molecule m;
    m.Insert(Site{16, Vector3D<double>(0,0,0)});
    m.Insert(Site{1, Vector3D<double>(0, r*std::sin(h), r*std::cos(h))});
    m.Insert(Site{1, Vector3D<double>(0,-r*std::sin(h), r*std::cos(h))});
    return m;
}
static void ExpectSADMatchesCoreGuess(const Molecule& mol, const std::string& basis, double tol)
{
    const auto sad=MoleculeSetup::Resolve(mol, {.model=Model::LDA, .molecule={.basis=basis}});
    ASSERT_EQ(sad->opts.seed, ChargeDensity::SeedStrategy::SAD) << "every element is in the all-electron library, so the Auto rule seeds from SAD";
    auto a=MoleculeCalculation::Prepare(sad).TakeValue(); (void)a->Converge();
    auto b=MoleculeCalculation::Prepare(MoleculeSetup::Resolve(mol, {.model=Model::LDA, .seed=ChargeDensity::SeedStrategy::CoreGuess, .molecule={.basis=basis}})).TakeValue(); (void)b->Converge();
    ASSERT_TRUE(a->Result()) << "the SAD start converges";
    ASSERT_TRUE(b->Result()) << "the core-guess start converges";
    EXPECT_NEAR(a->Energy(), b->Energy(), tol) << "one converged state from either start";
}
TEST(M_DFT, H2S_AllElectron_SADSeed_eqCoreGuess)    { ExpectSADMatchesCoreGuess(MakeH2S(),  "dzvp", 1e-6); }
TEST(M_DFT, TeF6_AllElectron_SADSeed_eqCoreGuess_Long) { ExpectSADMatchesCoreGuess(MakeTeF6(), "dzvp", 1e-6); }

TEST(M_DFT, WaterLDA_AuxiliaryFit)
{
    auto calc=MoleculeCalculation::Prepare(MoleculeSetup::Resolve(qchem::StructureData::GetMolecule("H2O"), {.model=Model::LDA, .vxcFit=Hamiltonian::VxcFit::Auxiliary, .coulomb={.auxBasis="A1_coul"}, .molecule={.basis="dzvp"}})).TakeValue(); (void)calc->Converge();
    ASSERT_TRUE(calc->Result());
    // The Hartree is STATED as A1_coul so the oracle (PySCF's A1_coul J-fit route) stays like-for-like: the δ default (user ruling 2026-11-22)
    // would move this cell by the A1 Coulomb fit's own 2.1 mHa.
    // RE-JUDGED 2026-11-15 (pin 10, lane S): FittedVxc retired.  The cell's E is now the exact XC of rho[D] on the run's mesh (FittedVxc fitted eps
    // for E too, the -1.4 mHa it sat off PySCF); only H still sees the A1 fit of v, a second-order error in E.  From -75.8808368818.  Independent
    // route: PySCF's exact XC, 2.1 μHa away (was 1.39 mHa).
    EXPECT_NEAR(calc->Energy(), -75.8794494584, 1e-5);   // anchor
    EXPECT_NEAR(calc->Energy(), -75.87944736, 2e-5) << "ORACLE: PySCF A1_coul J-fit route, exact XC; the gap is the A1 fit of v in H (second order)";
}

TEST(M_DFT, OxygenTripletLDA_AuxiliaryFit)
{
    auto calc=MoleculeCalculation::Prepare(MoleculeSetup::Resolve(MakeO2(), {.multiplicity=3, .model=Model::LDA, .vxcFit=Hamiltonian::VxcFit::Auxiliary, .coulomb={.auxBasis="A1_coul"}, .scf=kTight80, .molecule={.basis="dzvp"}})).TakeValue();
    ASSERT_TRUE(calc->Converge());
    // RE-JUDGED 2026-11-15 (pin 10, lane S): FittedVxc retired (see WaterLDA_AuxiliaryFit).  From -149.3018824495.  Independent routes: our δ cell
    // (OxygenTripletLDA) -149.3006843397, 7.5 μHa BELOW (the variational minimum: the fitted v only perturbs D), and PySCF 2.2 μHa away.
    EXPECT_NEAR(calc->Energy(), -149.3006768238, 1e-4);   // anchor
    EXPECT_NEAR(calc->Energy(), -149.30067460, 2e-5) << "ORACLE: PySCF UKS A1_coul J-fit route, exact XC; the gap is the A1 fit of v in H (second order)";
}

// F2 GTH q7 (the case the "Molecular PP route is 4 mHa/atom off CP2K" row was filed on) on BOTH routes: the δ default is CP2K's
// number to 1.4 mHa (PySCF's own A1_coul J-fit route to 57 μHa) -- the PP route was never the defect.
TEST(M_DFT, F2_PP)
{
    auto calc=MoleculeCalculation::Prepare(MoleculeSetup::Resolve(MakeF2(), {.model=Model::LDA, .pp="gth-lda-valence", .molecule={.basis="valence_lowq_vb"}})).TakeValue(); (void)calc->Converge();
    ASSERT_TRUE(calc->Result());
    // RE-JUDGED 2026-11-22 (pin 10, lane S; user ruling "Yes, δ default"): the Hartree is the δ fit.  From -48.2018945557 (the A1_coul Coulomb
    // fit, 1.35 mHa low); PySCF EXACT J and XC -48.20054460, 1.6 μHa away; CP2K -48.20058, 37 μHa.  (2026-11-05: the closed-form PP
    // integrals moved it +4.2e-5.)
    EXPECT_NEAR(calc->Energy(), -48.2005429519, 1e-5);   // anchor (δ default)
    EXPECT_NEAR(calc->Energy(), -48.20054460, 2e-5) << "ORACLE: PySCF, exact J and XC (CP2K -48.20058)";
}
TEST(M_DFT, F2_PP_AuxiliaryFit)
{
    auto calc=MoleculeCalculation::Prepare(MoleculeSetup::Resolve(MakeF2(), {.model=Model::LDA, .pp="gth-lda-valence", .vxcFit=Hamiltonian::VxcFit::Auxiliary, .coulomb={.auxBasis="A1_coul"}, .molecule={.basis="valence_lowq_vb"}})).TakeValue(); (void)calc->Converge();
    ASSERT_TRUE(calc->Result());
    // RE-JUDGED 2026-11-15 (pin 10, lane S): FittedVxc retired (see WaterLDA_AuxiliaryFit).  From -48.1970563871 (+4.8 mHa off PySCF: the eps fit
    // of a PP valence density).  Independent route: our δ cell (F2_PP) -48.2018945557, 93 μHa BELOW (variational); PySCF 78 μHa away.  The
    // larger second-order gap than water's is A1_exch's poorer fit of v on a pseudo-valence density (its P-rho is 29% negative, OpenWork).
    EXPECT_NEAR(calc->Energy(), -48.2018018561, 1e-5);   // anchor
    EXPECT_NEAR(calc->Energy(), -48.20187971, 1e-4) << "ORACLE: PySCF A1_coul J-fit route, exact XC; the gap is the A1 fit of v in H (second order)";
}

// THE A1 Pρ CELL IS RETIRED (user 2026-11-15, pin 28): XC must see ρ[D], and the A1 fit's projection of it is not ρ[D] (water's Pρ held 9.958
// of 10 electrons; F₂'s PP valence density went negative over 29% of the points).  Its last anchors (water -75.8990130799, F₂ PP
// -48.2067928423) and their PySCF first-order Pρ oracle (scripts/pyscf_prho_oracle.py) are in doc/OOD-SOLID-Cleanup.md V-XC-ONE-TERM.
// The setup now refuses the cell by name; the unit gate CalculationFacade.TheGaussianDeltaXCCellWarnsAndRefusesGDM holds the refusal.
TEST(M_DFT, F2_PP_AuxiliaryFitProjectedIsRefused)
{
    try
    {
        MoleculeSetup::Resolve(MakeF2(), {.model=Model::LDA, .pp="gth-lda-valence", .vxcFit=Hamiltonian::VxcFit::Auxiliary, .xcRho=Hamiltonian::XCρ::Fit, .molecule={.basis="valence_lowq_vb"}});
        ADD_FAILURE() << "xc.rho fit on the A1 fit ran";
    }
    catch (const std::invalid_argument& e) { EXPECT_NE(std::string(e.what()).find("pin 28"), std::string::npos) << e.what(); }
}

// THE δ HARTREE ROUTE (CoulombFit::Delta): ρ on the XC mesh, the direct Coulomb sum, first-order energy, through the same FittedVee.
// The direct sum is coarse, by measurement (TeF₆ prototype: J to about 1e-2), so the oracle tolerance says so.
TEST(M_DFT, WaterLDA_DeltaCoulomb_Long)
{
    auto calc=MoleculeCalculation::Prepare(MoleculeSetup::Resolve(qchem::StructureData::GetMolecule("H2O"), {.model=Model::LDA, .coulomb={.fit = Hamiltonian::CoulombFit::Delta, .solver = Hamiltonian::CoulombSolver::DirectSum}, .molecule={.basis="dzvp"}})).TakeValue(); (void)calc->Converge();
    ASSERT_TRUE(calc->Result());
    // The mesh pilot sizes XC and leaves this Hartree out (its error does not fall with the mesh, HartreeConvergesWithMesh): chasing it
    // walked the ladder at O(N²) per rung and ran ~10 h (2026-11-08).  A settled climb is the proof it was not chased.
    EXPECT_TRUE(calc->PilotCost().settled);
    EXPECT_NEAR(calc->Energy(), -75.8917526160, 1e-5);   // anchor
    EXPECT_NEAR(calc->Energy(), -75.87729931, 2e-2) << "ORACLE: PySCF, exact J and XC; the gap (-14.5 mHa) is the direct sum";
}

// ...and with the radial Poisson solve, the route to use: the same mesh, the operator at the mesh's own limit.  It sits 15x closer
// to exact J than the Gaussian A1 fit does (-2.1 mHa, the PySCF A1_coul row).  RE-JUDGED 2026-11-08 (pin 10, lane B): the unstated mesh
// now climbs to eps (GeometryOptimizationPlan §14): δ-ρ Hartree rides it and misses eps x N at (80, 29), so the pilot raises it to
// (100, 29), and the energy moves -39.9 μHa, from +30.6 to -9.3 μHa off the exact-J oracle.
TEST(M_DFT, WaterLDA_DeltaCoulomb_Poisson)
{
    auto calc=MoleculeCalculation::Prepare(MoleculeSetup::Resolve(qchem::StructureData::GetMolecule("H2O"), {.model=Model::LDA, .coulomb={.fit = Hamiltonian::CoulombFit::Delta}, .molecule={.basis="dzvp"}})).TakeValue(); (void)calc->Converge();
    ASSERT_TRUE(calc->Result());
    EXPECT_NEAR(calc->Energy(), -75.8773086479, 1e-6);   // anchor (re-judged 2026-11-08, above)
    EXPECT_NEAR(calc->Energy(), -75.87729931, 5e-5) << "ORACLE: PySCF, exact J and XC; the gap (-9 μHa) is the mesh's, inside eps x N";
}

// The multipole fit has δ x Poisson's potential, so its SCF is the same one; only the energy's second leg differs.  Re-judged with the
// δ case above (2026-11-08): the pilot's (100, 29) moves it -39.5 μHa, from +32 to -7.9 μHa off the oracle.
TEST(M_DFT, WaterLDA_MultipoleCoulomb)
{
    auto calc=MoleculeCalculation::Prepare(MoleculeSetup::Resolve(qchem::StructureData::GetMolecule("H2O"), {.model=Model::LDA, .coulomb={.fit = Hamiltonian::CoulombFit::Multipole}, .molecule={.basis="dzvp"}})).TakeValue(); (void)calc->Converge();
    ASSERT_TRUE(calc->Result());
    EXPECT_NEAR(calc->Energy(), -75.8773072275, 1e-6);   // anchor (re-judged 2026-11-08, above)
    EXPECT_NEAR(calc->Energy(), -75.87729931, 5e-5) << "ORACLE: PySCF, exact J and XC: delta's level, -8 μHa.  This rung was 530 out when rho-tilde was synthesized onto the mesh";
}

// The multipole route's energy depends less on the l cut than δ's: at 150/35, 11 μHa between l = 12 and 16 where δ moves 35.
TEST(M_DFT, WaterLDA_MultipoleCoulomb_eqLCut_Long)
{
    auto energy=[](int lMax)
    {
        auto calc=MoleculeCalculation::Prepare(MoleculeSetup::Resolve(qchem::StructureData::GetMolecule("H2O"),
            {.model=Model::LDA, .xcMesh=qcMesh::MeshParams{.radial = qcMesh::RadialKind::MHL, .nRadial = 150, .mhl_m = 3, .mhl_alpha = 2.0,
                                                      .angular = qcMesh::AngularKind::Lebedev, .angularDegree = 35, .beckeOrder = 2}, .coulomb={.fit = Hamiltonian::CoulombFit::Multipole, .lMax = lMax}, .molecule={.basis="dzvp"}})).TakeValue(); (void)calc->Converge();
        EXPECT_TRUE(calc->Result());
        return calc->Energy();
    };
    const double e12=energy(12);
    EXPECT_NEAR(e12, energy(16), 2e-5);
    EXPECT_NEAR(e12, -75.87729931, 2e-5) << "ORACLE: PySCF, exact J and XC";
}

// Multipole for Vxc, through the fitter as the A1 route is.  The Hartree stays on the default route, so only the XC term changed.
// RE-JUDGED 2026-11-22 (pin 10, lane S; user ruling "Yes, δ default"): the default Hartree is the δ fit now, so the oracle is PySCF's EXACT J
// (-75.87729931); this run sits -9.3 μHa off it (it sat -2.14 mHa off it, i.e. on the A1_coul J-fit value, before the flip).
// RE-JUDGED 2026-11-08 (pin 10, lane B): the gap was the MESH's, not the multipole Vxc's.  At (80, 29) this run sat -140 μHa off the
// A1-Hartree value; the mesh pilot (GeometryOptimizationPlan §14) climbs it to (100, 41 GL), and it sits -5.1 μHa off, inside eps x N.
TEST(M_DFT, WaterLDA_MultipoleVxc)
{
    auto calc=MoleculeCalculation::Prepare(MoleculeSetup::Resolve(qchem::StructureData::GetMolecule("H2O"), {.model=Model::LDA, .vxcFit=Hamiltonian::VxcFit::Multipole, .molecule={.basis="dzvp"}})).TakeValue(); (void)calc->Converge();
    ASSERT_TRUE(calc->Result());
    EXPECT_NEAR(calc->Energy(), -75.87729931, 3e-5) << "ORACLE: PySCF, exact J and XC; eps x N_atoms, the run's own target";
}

// --- Ionic SAD on a molecule (user, 2026-11-12: "we should strive to make seeding ... structure neutral") -------------------------
// The SAD superposition at each site's formal charge, from the ONE charge rule the periodic seed uses (ChargeDensity::IonicFormalCharges).
// A seed is a starting point, so the claim is that it reaches the SAME state as the neutral SAD seed: a closed-shell covalent molecule
// (water: the rule makes it H+ O2-, the hardest case for a seed), an ionic one (LiF: Li+ F-), and a charged one (OH-, whose extra
// electron only the fit's charge constraint carries).  The iteration counts print, so a seed that converges slower is visible.
namespace
{
Molecule MakeLiF()    // Li-F 1.564 Angstrom = 2.9555 bohr, along z
{
    Molecule m;
    m.Insert(Site{3, Vector3D<double>(0, 0, 0)});
    m.Insert(Site{9, Vector3D<double>(0, 0, 2.9555)});
    return m;
}
Molecule MakeHydroxide()    // OH-: O-H 0.97 Angstrom = 1.8330 bohr; the anion's extra electron on the O site
{
    Molecule m;
    m.Insert(Site{8, Vector3D<double>(0, 0, 0), -1.0});
    m.Insert(Site{1, Vector3D<double>(0, 0, 1.8330)});
    return m;
}
void IonicSADReachesTheSADState(const Molecule& mol, const char* what)
{
    using qchem::ChargeDensity::SeedStrategy;
    double E[2]={0,0};
    const SeedStrategy seeds[2]={SeedStrategy::SAD, SeedStrategy::IonicSAD};
    for (int i=0;i<2;i++)
    {
        auto calc=MoleculeCalculation::Prepare(MoleculeSetup::Resolve(mol, {.model=Model::LDA, .seed=seeds[i], .molecule={.basis="dzvp"}})).TakeValue(); (void)calc->Converge();
        const auto r=calc->Result();
        ASSERT_TRUE(r) << what << (i ? " IonicSAD" : " SAD") << " did not converge";
        E[i]=r->Energy();
        std::cout << "[ionic seed] " << what << (i ? " IonicSAD: " : " SAD:       ") << r->IterationCount() << " iterations, E = "
                  << std::setprecision(12) << E[i] << std::endl;
    }
    EXPECT_NEAR(E[1], E[0], 1e-7) << what << ": the ionic seed must reach the SAD state";
}
}
TEST(M_DFT, WaterLDA_IonicSADReachesTheSADState)     { IonicSADReachesTheSADState(qchem::StructureData::GetMolecule("H2O"), "H2O"); }
TEST(M_DFT, LiF_LDA_IonicSADReachesTheSADState)      { IonicSADReachesTheSADState(MakeLiF(), "LiF"); }
TEST(M_DFT, HydroxideLDA_IonicSADReachesTheSADState) { IonicSADReachesTheSADState(MakeHydroxide(), "OH-"); }

// A CHARGED molecule under a PSEUDOPOTENTIAL keeps its charge (lane G, 2026-11-21): MoleculeSetup's point structure REPLACED each site's
// stated net charge with Z - Zion, so OH- with GTH ran as NEUTRAL OH (7 valence electrons), silently.  The oracle is PySCF 2.14 on the same
// model: CP2K's DZVP-MOLOPT-GTH read from BASIS_MOLOPT (spherical), GTH-PADE, LDA "lda,vwn" (Dirac + VWN5), grids level 7, conv 1e-11:
// OH- (8 e) -16.5335495587; the uncharged run (7 e) landed at -16.4651900 on today's code before the fix.  The Hartree term runs the
// DELTA route (rho on the mesh, 150/47) because the default A1 Coulomb fit is too small for this basis, 1.9 mHa LOW (4.4 on H2O; OpenWork
// row "The Coulomb-fit default is a hidden policy"); on the delta route ours is 0.4 uHa from PySCF (H2O at the same settings: 1.3 uHa).
TEST(M_DFT, Hydroxide_PP_KeepsItsChargeVsPySCF)
{
    RunOptions o{.model=Model::LDA, .pp="gth-lda-valence", .molecule={.basis="dzvp_molopt_gth"}};
    o.molecule.angular=Angular::Spherical;
    o.coulomb.fit=Hamiltonian::CoulombFit::Delta;
    MoleculeXCMesh(o).nRadial=150; MoleculeXCMesh(o).angularDegree=47;
    auto calc=MoleculeCalculation::Prepare(MoleculeSetup::Resolve(MakeHydroxide(), o)).TakeValue(); (void)calc->Converge();
    ASSERT_TRUE(calc->Result()) << "OH- with GTH must converge";
    std::cout << "[charged PP] OH- E = " << std::setprecision(12) << calc->Energy() << std::endl;
    EXPECT_NEAR(calc->Energy(), -16.5335495587, 1e-5) << "PySCF, same basis / PP / functional";
}
