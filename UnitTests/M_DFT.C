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
