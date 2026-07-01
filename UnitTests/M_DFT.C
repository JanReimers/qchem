// File: UnitTests/M_DFT.C  Molecular un-polarized Xalpha-DFT total-energy tests (facade-driven).
//
// Production path: qchem::Calculation builds the BasisSet::Molecule::Factory PolarizedGaussian (M&D)
// basis; the Xalpha DFT Hamiltonian adds the fitted Coulomb + Xalpha exchange on a numerical mesh.  We
// drive a full SCF for N2 and water through the facade and check the converged total energy against a
// reference -- the same total-energy check used for the atom A_DFT tests.
//
// Migrated off the QchemTester/TestMolecule scaffold (OpenWork D).  The facade auto-selects the SAD seed
// AND DIIS-from-start for DFT (CalcOptions: model is DFT -> EMax=100 so the [F,D] gate never blocks DIIS),
// which is exactly the DIIS_FromStart + SAD recipe the scaffold set by hand here -- so anchors are
// byte-for-byte the scaffold values; the default mesh matches TestMolecule::GetMeshParams() exactly.
#include "gtest/gtest.h"
#include <cmath>
#include <stdexcept>

import qchem.Calculation;            // Calculation, CalcOptions, Model, Angular
import qchem.Structure;              // Molecule, Atom
import qchem.Types;                  // Vector3D
import qchem.Unittests.TestUtils;    // RelativeError
using namespace qchem;

static Molecule MakeN2()
{
    Molecule m;
    m.Insert(new Atom(7, 0, Vector3D<double>(-1.03, 0, 0)));
    m.Insert(new Atom(7, 0, Vector3D<double>( 1.04, 0, 0)));
    return m;
}
static Molecule MakeWater()      // experimental geometry in BOHR, C2 axis along z
{
    Molecule w;
    w.Insert(new Atom(8, 0, Vector3D<double>(0,  0.0,   0.0)));
    w.Insert(new Atom(1, 0, Vector3D<double>(0,  1.431, 1.107)));
    w.Insert(new Atom(1, 0, Vector3D<double>(0, -1.431, 1.107)));
    return w;
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
    Calculation calc(MakeN2(), {.basis = "dzvp", .model = Model::Xalpha, .xalpha = 0.75197});
    EXPECT_LT(fabs(RelativeError(calc.Energy(), -109.21679)), 2e-3);
}

TEST(M_DFT, Water)
{
    Calculation calc(MakeWater(), {.basis = "dzvp", .model = Model::Xalpha, .xalpha = 0.74000});
    // CONVERGED Xalpha total energy (the facade's auto DIIS-from-start): [F,D] and Δρ reach ~1e-13 by
    // ~iter 14, virial V/K = -2.006 (the correct ~-2, vs the buggy -2.11).  A true converged regression
    // sentinel.  (Was -76.123348 interim limit-cycle, and -79.414120 stale-cache.)
    EXPECT_LT(fabs(RelativeError(calc.Energy(), -76.1493013984)), 2e-3);
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
    Calculation calc(MakeWater(), {.basis   = "dzvp", .model   = Model::Xalpha,
                                   .xalpha  = 0.74000, .angular = Angular::Spherical});
    EXPECT_NEAR(calc.Energy(), -76.1493013984, 0.05);                 // sanity: same ballpark as Cartesian
    EXPECT_LT(fabs(RelativeError(calc.Energy(), -76.1485556624)), 2e-3);  // CONVERGED anchor (DIIS from start)
}

// Polarized Xalpha water through the facade's DEFAULT (auto SAD) seed.  This used to SEGFAULT: the SAD seed
// is a spin-agnostic total density, and the polarized Vxc null-derefed when its dynamic_cast<Polarized_CD>
// of the seed failed (the assert was compiled out in Release).  Fixed in FittedVxcPol::CalcMatrix -- a
// spin-unpolarized seed (rho_up=rho_down=rho/2) maps each channel's iteration-0 Vxc to the unpolarized Vxc
// of the total density.  Water is closed shell, so with the SAME alpha the polarized run must converge to
// the SAME energy as the unpolarized anchor -- and it does, to ~1e-11 (confirming a correct seed Fock).
TEST(M_DFT, WaterPolarizedSAD)
{
    Calculation calc(MakeWater(), {.basis = "dzvp", .model = Model::Xalpha, .xalpha = 0.74000, .pol = Pol::Polarized});
    // Polarized water Xalpha is the oscillatory case (commutator ~2-4; see M_Sym RunDFT) -- it needs more
    // iterations + damping than the facade's quick default, so re-Converge with tight params (DIIS already
    // drives from the start for DFT).
    calc.Converge({.NMaxIter = 60, .MinΔρ = 1e-7, .MinΔFD = 1e-9, .MinVirial = 1e2, .MinFD = 1e-7,
                   .StartingRelaxRo = 0.5, .MergeTol = 1e-4, .Verbose = false});
    EXPECT_LT(fabs(RelativeError(calc.Energy(), -76.1493013984)), 2e-3);   // == the unpolarized Xalpha anchor
}

// Spin-native (polarized) parameter-free LDA water (Dirac exchange + spin-native VWN5 correlation) through
// the facade's auto SAD seed -- OpenWork B2: Ham_DFTcorr_P + FittedVcorrPol.  Water is closed shell, so at
// convergence zeta=0 everywhere and the polarized run must collapse to the UNPOLARIZED LDA anchor
// (M_Calculation.WaterLDA = -75.9324615507).  This exercises the new polarized correlation term end to end:
// the SAD-seed fallback (rho_up=rho_down=rho/2), the per-spin v_c^sigma fit, and the two-channel E_c.
TEST(M_DFT, WaterPolarizedLDA)
{
    Calculation calc(MakeWater(), {.basis = "dzvp", .model = Model::LDA, .pol = Pol::Polarized});
    calc.Converge({.NMaxIter = 60, .MinΔρ = 1e-7, .MinΔFD = 1e-9, .MinVirial = 1e2, .MinFD = 1e-7,
                   .StartingRelaxRo = 0.5, .MergeTol = 1e-4, .Verbose = false});
    EXPECT_NEAR(calc.Energy(), -75.9324615507, 1e-6);   // collapses to the unpolarized LDA anchor (closed shell)
}

// OPEN-SHELL through the facade multiplicity knob (OpenWork B4): triplet O2 (16 e-, multiplicity 3 =>
// nUp=9, nDown=7).  multiplicity>1 auto-promotes the calc to Pol::Polarized, builds Molecule_EC(9,7), and
// runs spin-native LDA.  The spin texture is real here (the half-filled pi* gives the triplet), so this
// genuinely exercises zeta!=0 -- not the closed-shell collapse the WaterPolarized tests hit.  "Did E move"
// regression sentinel (converged spin-native LDA total energy).
TEST(M_DFT, OxygenTripletLDA)
{
    Calculation calc(MakeO2(), {.basis = "dzvp", .model = Model::LDA, .multiplicity = 3});
    calc.Converge({.NMaxIter = 80, .MinΔρ = 1e-7, .MinΔFD = 1e-9, .MinVirial = 1e2, .MinFD = 1e-7,
                   .StartingRelaxRo = 0.5, .MergeTol = 1e-4, .Verbose = false});
    EXPECT_NEAR(calc.Energy(), -149.2562876393, 1e-4);   // converged spin-native LDA triplet (regression anchor)
}

// Physics check that the spin channel is doing real work: the triplet ground state must lie BELOW the
// closed-shell singlet (Hund's rule).  A flat/zero gap would mean the open-shell occupation silently
// collapsed to closed shell -- the failure mode B3/B4 exist to prevent.
TEST(M_DFT, OxygenTripletBelowSinglet)
{
    auto run = [](int mult) {
        Calculation c(MakeO2(), {.basis = "dzvp", .model = Model::LDA, .multiplicity = mult});
        c.Converge({.NMaxIter = 80, .MinΔρ = 1e-7, .MinΔFD = 1e-9, .MinVirial = 1e2, .MinFD = 1e-7,
                    .StartingRelaxRo = 0.5, .MergeTol = 1e-4, .Verbose = false});
        return c.Energy();
    };
    const double triplet = run(3);   // nUp=9, nDown=7  (Pol::Polarized, auto)
    const double singlet = run(1);   // nUp=8, nDown=8  (closed shell)
    EXPECT_LT(triplet, singlet - 0.01);   // triplet is the ground state, by a clear margin
}

// A multiplicity whose parity disagrees with Ne must be rejected, not silently miscounted: water has 10
// electrons (even), so a doublet (multiplicity 2 => 2S=1, odd) is impossible.  Fail loud (B4 validation).
TEST(M_DFT, BadMultiplicityThrows)
{
    EXPECT_THROW(Calculation(MakeWater(), {.basis = "sto-3g", .model = Model::LDA, .multiplicity = 2}),
                 std::runtime_error);
}
