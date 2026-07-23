// File: UnitTests/M_HF_U.C  Molecular un-polarized Hartree-Fock total-energy tests (facade-driven).
//
// Production path: qchem::Calculation builds the BasisSet::Molecule::Factory PolarizedGaussian (M&D)
// basis and runs the canonical SCF.  Here we drive a full SCF for N2 and water through the facade and
// check the converged total energy against a reference -- the same total-energy check used for the atom
// A_HF tests.  (Integral correctness to machine precision is covered separately by M_PG_Oracle.)
//
// Migrated off the QchemTester/TestMolecule scaffold (OpenWork D): the facade IS the production recipe,
// so the tests now read like end-user code.  Anchors are byte-for-byte the scaffold values -- the
// facade's default SCFParams are identical to the old `scf` literal, so convergence is unchanged.
#include "gtest/gtest.h"
#include <cmath>

import qchem.Calculation;            // Calculation, CalcOptions (+ Model/Pol/Engine/Angular)
import qchem.Structure;              // Molecule, Atom
import qchem.Types;                  // Vector3D
import qchem.ChargeDensity.Seed;     // SeedStrategy
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

#ifndef DEBUG
// Anchors are the converged HF/DZVP total energies (pinned to catch movement -- "did E move",
// not a physical-accuracy claim).
TEST(M_HF_U, N2)            // slow in Debug (~80s); Release-only
{
    Calculation calc(MakeN2(), {.basis = "dzvp"});
    EXPECT_LT(fabs(RelativeError(calc.Energy(), -108.94524)), 1e-4);
}
#endif

TEST(M_HF_U, Water)
{
    Calculation calc(MakeWater(), {.basis = "dzvp"});
    EXPECT_LT(fabs(RelativeError(calc.Energy(), -76.022903)), 1e-5);
}

// --- HF from a SAD seed: the matrix-free-seed bootstrap --------------------------------------------
// HF's exact-exchange K needs the density MATRIX, so a matrix-free SAD (superposition-of-atomic-densities)
// seed can't build the iteration-0 Fock directly.  The SCFIterator bootstraps it: runs the SAD rho through
// a default LDA DFT sibling for one step to manufacture a real D0, then seeds HF with D0.  The converged
// energy is seed-independent, so it must land on the SAME CoreGuess anchor -- that equality is the proof
// the bootstrap produced a valid density.  See project_numericcd_refactor.
TEST(M_HF_U, WaterSADseed)
{
    Calculation calc(MakeWater(), {.basis = "dzvp", .seed = ChargeDensity::SeedStrategy::SAD});
    EXPECT_LT(fabs(RelativeError(calc.Energy(), -76.022903)), 1e-5);  // identical converged E to CoreGuess
}

// --- The Cartesian basis integrated by libcint through the same SCF --------------------------------
// engine=LibCint selects the PG_LibCint matrix-delivery evaluator (same dzvp data, same Cartesian
// component set as the M&D path).  Since libcint's integrals match M&D to <1e-10 element-wise (M_LibCint
// guard), the converged HF energy must equal the Cartesian anchor -- the end-to-end cross-check of the
// external engine through a real SCF.
TEST(M_HF_U, WaterLibCint)
{
    Calculation calc(MakeWater(), {.basis = "dzvp", .engine = Engine::LibCint});
    EXPECT_LT(fabs(RelativeError(calc.Energy(), -76.022903)), 1e-5);  // same converged E as the M&D path
}

// libcint with native real-spherical functions (engine=libcint + spherical) -- the INDEPENDENT spherical
// oracle PG_Spherical lacked.  The HF energy is basis-ordering invariant, so even though libcint's harmonic
// order/convention differs from PG_Spherical's it must converge to the same -76.020277 (the M&D spherical
// value), validating the spherical basis end-to-end against a second engine.
TEST(M_HF_U, WaterLibCintSpherical)
{
    Calculation calc(MakeWater(), {.basis   = "dzvp",
                                   .engine  = Engine::LibCint,
                                   .angular = Angular::Spherical});
    EXPECT_LT(fabs(RelativeError(calc.Energy(), -76.020277)), 1e-5);  // == PG_Spherical (M&D) spherical HF
}

// --- The spherical-Gaussian (PG_Spherical) basis through the same SCF ------------------------------
// angular=Spherical selects the real-solid-harmonic basis.  The 5 spherical d's are a strict subset of
// the 6 Cartesian d's (the missing 6th, xx+yy+zz, is an l=0 contaminant), so the spherical basis is
// variationally a SUBSET: E_sph >= E_cart, by only the small contaminant gap.  That inequality
// (oracle-free) validates the new basis end-to-end through the production SCF; the pinned value is the
// usual "did E move" regression anchor.
TEST(M_HF_U, WaterSpherical)
{
    Calculation calc(MakeWater(), {.basis = "dzvp", .angular = Angular::Spherical});
    const double E_cart = -76.022903;          // the Cartesian-DZVP anchor (water has one O d shell)
    const double E_sph  = calc.Energy();
    EXPECT_GT(E_sph, E_cart);                  // variational subset: dropping the d contaminant raises E
    EXPECT_LT(E_sph - E_cart, 0.05);           // ... but only slightly (the contaminant is a minor fn)
    EXPECT_LT(fabs(RelativeError(E_sph, -76.020277)), 1e-5);  // regression anchor (just above E_cart)
}
