// File: UnitTests/M_DFT.C  Molecular un-polarized Xalpha-DFT total-energy tests (facade-driven).
//
// Production path: qchem::Calculation builds the BasisSet::Gaussian::Factory PolarizedGaussian (M&D)
// basis; the Xalpha DFT Hamiltonian adds the fitted Coulomb + Xalpha exchange on a numerical mesh.  We
// drive a full SCF for N2 and water through the facade and check the converged total energy against a
// reference -- the same total-energy check used for the atom A_DFT tests.
//
// Migrated off the QchemTester/TestMolecule scaffold (OpenWork D).  The facade auto-selects the SAD seed
// AND DIIS-from-start for DFT (CalcOptions: model is DFT -> EMax=100 so the [F,D] gate never blocks DIIS),
// which is exactly the DIIS_FromStart + SAD recipe the scaffold set by hand here -- so anchors are
// byte-for-byte the scaffold values.
//
// RE-ANCHORED 2026-10-31 (lane B) on the CONVERGED default XC mesh (MHL 80 / Lebedev 29; was 30 / 5, 20 mHa/atom low on
// water -- MoleculeSetup.C).  Re-judged per pin 10 against PySCF (exact J and XC, same Cartesian dzvp, 200/590 grid):
// every closed-shell anchor now sits 1-5 mHa ABOVE it -- the DGauss A1 fitted-Coulomb/Vxc route (PySCF's A1_coul J-fit alone
// moves 2 mHa DOWN, so the fitted Vxc carries it) -- where the 30/5 values sat 10-55 mHa BELOW (not variational).
#include "gtest/gtest.h"
#include <cmath>
#include <stdexcept>

import qchem.Calculation;            // Calculation, CalcOptions, Model, Angular
import qchem.SCFCap;                             // IterationBudget: the convergence PANEL's budget (doc/SCFRobustnessPlan.md §4)
import qchem.SCFParams;              // SCFParams (the O2 runs converge once, on their own recipe)
import qchem.StructureData;        // GetMolecule("H2O") -- the shared water geometry (D-MAKEWATER)
import qchem.Structure;              // Molecule, Atom
import qchem.Types;                  // Vector3D
import qchem.PeriodicTable;          // RelativeError
using namespace qchem;

static Molecule MakeN2()
{
    Molecule m;
    m.Insert(new Atom(7, 0, Vector3D<double>(-1.03, 0, 0)));
    m.Insert(new Atom(7, 0, Vector3D<double>( 1.04, 0, 0)));
    return m;
}
static Molecule MakeO2()         // O2 at the experimental bond length ~1.208 Ang = 2.282 bohr, along z
{
    Molecule m;
    m.Insert(new Atom(8, 0, Vector3D<double>(0, 0,  1.141)));
    m.Insert(new Atom(8, 0, Vector3D<double>(0, 0, -1.141)));
    return m;
}

// Anchors are the converged Xalpha total energies (deterministic regression sentinels -- "did E
// move", NOT the physical energy; exchange-only Xalpha at these alphas is not calibrated to it).
TEST(M_DFT, N2)
{
    Calculation calc(MoleculeSetup::Resolve(MakeN2(), {.basis = "dzvp", .model = Model::Xalpha, .xalpha = 0.75197}));
    // -109.225960 at the converged mesh (MHL 150/35: -109.225954); PySCF exact -109.227199 (+1.2 mHa, the fitted Vxc).
    // Was -109.21679 at 2e-3 relative (0.2 Ha!) on the 30/5 mesh (-109.217602): the tolerance hid the mesh.
    EXPECT_NEAR(calc.Energy(), -109.2259603700, 1e-5);
}

TEST(M_DFT, Water)
{
    Calculation calc(MoleculeSetup::Resolve(qchem::StructureData::GetMolecule("H2O"), {.basis = "dzvp", .model = Model::Xalpha, .xalpha = 0.74000}));
    // CONVERGED Xalpha total energy (the facade's auto DIIS-from-start): [F,D] and Δρ reach ~1e-13 by
    // ~iter 14, virial V/K = -2.006 (the correct ~-2, vs the buggy -2.11).  A true converged regression
    // sentinel.  (Was -76.123348 interim limit-cycle, and -79.414120 stale-cache, and -76.149301 on the 30/5 mesh -- 48 mHa low,
    // inside the old 2e-3 relative.)  Converged mesh (MHL 150/35) -76.100894; PySCF exact -76.104768 (+3.9 mHa, the fitted Vxc).
    EXPECT_NEAR(calc.Energy(), -76.1008971700, 1e-5);
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
    Calculation calc(MoleculeSetup::Resolve(qchem::StructureData::GetMolecule("H2O"), {.basis   = "dzvp", .model   = Model::Xalpha,
                                   .angular = Angular::Spherical, .xalpha  = 0.74000}));
    EXPECT_NEAR(calc.Energy(), -76.1008971700, 0.05);                 // sanity: same ballpark as Cartesian
    // CONVERGED anchor (DIIS from start).  MHL 150/35: -76.100117; PySCF exact (spherical dzvp) -76.104104 (+4.0 mHa, as Cartesian).
    // On the converged mesh the spherical lands 0.8 mHa ABOVE the Cartesian, not the 0.09 Ha the note above was written against.
    EXPECT_NEAR(calc.Energy(), -76.1001208265, 1e-5);
}

// Polarized Xalpha water through the facade's DEFAULT (auto SAD) seed.  This used to SEGFAULT: the SAD seed
// is a spin-agnostic total density, and the polarized Vxc null-derefed when its dynamic_cast<Polarized_CD>
// of the seed failed (the assert was compiled out in Release).  Fixed in FittedVxcPol::GetMatrix -- a
// spin-unpolarized seed (rho_up=rho_down=rho/2) maps each channel's iteration-0 Vxc to the unpolarized Vxc
// of the total density.  Water is closed shell, so with the SAME alpha the polarized run must converge to
// the SAME energy as the unpolarized anchor -- and it does, to ~1e-11 (confirming a correct seed Fock).
TEST(M_DFT, WaterPolarizedSAD)
{
    Calculation calc(MoleculeSetup::Resolve(qchem::StructureData::GetMolecule("H2O"), {.basis = "dzvp", .model = Model::Xalpha, .spin = SpinGroup::Polarized, .xalpha = 0.74000}));
    // Polarized water Xalpha is the oscillatory case (commutator ~2-4; see M_Sym RunDFT) -- it needs more
    // iterations + damping than the facade's quick default, so re-Converge with tight params (DIIS already
    // drives from the start for DFT).
    ASSERT_TRUE(calc.Converge({.NMaxIter = 60, .MinΔρPerElectron = 1e-7, .MinΔFD = 1e-9, .MinVirial = 1e2, .MinFD = 1e-7,
                   .StartingRelaxRo = 0.5, .MergeTol = 1e-4, .Verbose = false}));
    EXPECT_NEAR(calc.Energy(), -76.1008971700, 1e-5);   // == the unpolarized Xalpha anchor
}

// Spin-native (polarized) parameter-free LDA water (Dirac exchange + spin-native VWN5 correlation) through
// the facade's auto SAD seed -- OpenWork B2: Ham_DFTcorr_P + FittedVcorrPol.  Water is closed shell, so at
// convergence zeta=0 everywhere and the polarized run must collapse to the UNPOLARIZED LDA anchor
// (M_Calculation.WaterLDA = -75.8726517).  This exercises the new polarized correlation term end to end:
// the SAD-seed fallback (rho_up=rho_down=rho/2), the per-spin v_c^sigma fit, and the two-channel E_c.
TEST(M_DFT, WaterPolarizedLDA)
{
    Calculation calc(MoleculeSetup::Resolve(qchem::StructureData::GetMolecule("H2O"), {.basis = "dzvp", .model = Model::LDA, .spin = SpinGroup::Polarized}));
    ASSERT_TRUE(calc.Converge({.NMaxIter = 60, .MinΔρPerElectron = 1e-7, .MinΔFD = 1e-9, .MinVirial = 1e2, .MinFD = 1e-7,
                   .StartingRelaxRo = 0.5, .MergeTol = 1e-4, .Verbose = false}));
    EXPECT_NEAR(calc.Energy(), -75.8726517100, 1e-6);   // collapses to the unpolarized LDA anchor (closed shell)
}

// OPEN-SHELL through the facade multiplicity knob (OpenWork B4): triplet O2 (16 e-, multiplicity 3 =>
// nUp=9, nDown=7).  multiplicity>1 auto-promotes the calc to SpinGroup::Polarized, builds Molecule_EC(9,7), and
// runs spin-native LDA.  The spin texture is real here (the half-filled pi* gives the triplet), so this
// genuinely exercises zeta!=0 -- not the closed-shell collapse the WaterPolarized tests hit.  "Did E move"
// regression sentinel (converged spin-native LDA total energy).
TEST(M_DFT, OxygenTripletLDA)
{
    Calculation calc(MoleculeSetup::Resolve(MakeO2(), {.basis = "dzvp", .model = Model::LDA, .multiplicity = 3}));
    ASSERT_TRUE(calc.Converge({.NMaxIter = 80, .MinΔρPerElectron = 1e-7, .MinΔFD = 1e-9, .MinVirial = 1e2, .MinFD = 1e-7,
                   .StartingRelaxRo = 0.5, .MergeTol = 1e-4, .Verbose = false}));
    // Converged spin-native LDA triplet (regression anchor).  MHL 150/35: -149.263710; was -149.256288 on the 30/5 mesh.
    // ⚠ PySCF UKS lda,vwn (exact J and XC) gives -149.299178: we are 35.5 mHa ABOVE -- 8x water's fit offset, while our UHF
    // triplet matches PySCF to 1e-8.  An open row (OpenWork §3, "molecular polarized DFT is 35 mHa above PySCF on O₂").
    EXPECT_NEAR(calc.Energy(), -149.2637158314, 1e-4);
}

// REGRESSION for the O2-HF-triplet SEGV: the polarized level display (tCompositeWF::DisplayEigenPolarized) crashed on
// an open-shell level present in only ONE spin channel (find() is UB-in-Release on a miss).  Fixed via
// EnergyLevels::FindOrNull + labelling from the always-valid combined el.qns.  This runs it with Verbose ON
// (the crash path) and pins the converged HF triplet energy.  "Did E move" + no-crash.
TEST(M_DFT, OxygenTripletHF_VerboseDisplay)
{
    Calculation calc(MoleculeSetup::Resolve(MakeO2(), {.basis = "dzvp", .model = Model::HF, .multiplicity = 3}));
    ASSERT_TRUE(calc.Converge({.NMaxIter = 80, .MinΔρPerElectron = 1e-6, .MinΔFD = 1e-8, .MinVirial = 1e2, .MinFD = 1e-6,
                   .StartingRelaxRo = 0.5, .MergeTol = 1e-4, .Verbose = true}));   // Verbose => DisplayEigen (the crash path)
    EXPECT_NEAR(calc.Energy(), -149.6295470276, 1e-4);   // converged HF triplet O2
}

// Physics check that the spin channel is doing real work: the triplet ground state must lie BELOW the
// closed-shell singlet (Hund's rule).  A flat/zero gap would mean the open-shell occupation silently
// collapsed to closed shell -- the failure mode B3/B4 exist to prevent.
TEST(M_DFT, OxygenTripletBelowSinglet_Panel)   // PANEL case 11: a degenerate partly filled level (smeared singlet)
{
    qchem::IterationBudget budget(33, "O2 triplet + smeared singlet (molecular, degenerate pi*; 26 measured, 2026-10-29 panel baseline)");
    // The closed-shell SINGLET has two electrons for the exactly degenerate pi*x / pi*y pair: a single integer determinant
    // flips them between the two every iteration and never converged (2026-10-29, the cap guard).  A tiny smearing occupies
    // the pair EVENLY (equal fractions at any kT for a degenerate level) -- the symmetric closed-shell singlet, still above
    // the triplet by the exchange the triplet gains.  Integer occupations for the triplet (pi*x up, pi*y up: no tie).
    auto run = [](int mult) {
        SCFParams p{.NMaxIter = 100, .MinΔρPerElectron = 1e-7, .MinVirial = 1e2, .MinFD = 1e-7,
                    .StartingRelaxRo = 0.5, .MergeTol = 1e-4, .Verbose = false};
        if (mult==1) p.SmearingkT=1e-4;
        Calculation c(MoleculeSetup::Resolve(MakeO2(), {.basis = "dzvp", .model = Model::LDA, .multiplicity = mult}), p);
        EXPECT_TRUE(c.IsConverged()) << "multiplicity " << mult;
        return c.Energy();
    };
    const double triplet = run(3);   // nUp=9, nDown=7  (SpinGroup::Polarized, auto)
    const double singlet = run(1);   // nUp=8, nDown=8  (closed shell)
    EXPECT_LT(triplet, singlet - 0.01);   // triplet is the ground state, by a clear margin
}

// A multiplicity whose parity disagrees with Ne must be rejected, not silently miscounted: water has 10
// electrons (even), so a doublet (multiplicity 2 => 2S=1, odd) is impossible.  Fail loud (B4 validation).
TEST(M_DFT, BadMultiplicityThrows)
{
    EXPECT_THROW(Calculation(MoleculeSetup::Resolve(qchem::StructureData::GetMolecule("H2O"), {.basis = "sto-3g", .model = Model::LDA, .multiplicity = 2})),
                 std::runtime_error);
}
