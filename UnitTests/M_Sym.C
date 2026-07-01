// File: UnitTests/M_Sym.C  Molecular symmetry-adapted SCF integration tests (SALC end-to-end, facade-driven).
//
// Drives an SCF for water two ways through qchem::Calculation -- {.symmetry=false} (a plain single-IBS
// basis) and {.symmetry=true} (a SymmetryAdaptedBasisSet: one IBS per C2v irrep, global aufbau) -- and
// checks the total energies agree.  {.symmetry=true} runs exactly the production SALC path the scaffold
// used to wire by hand (PG::SymmetryAdapt over the raw PolarizedGaussian basis, Omega_ab shared across
// irreps via the global Cache2).  Covers HF and Xalpha-DFT (un/polarized) plus rigid rotation/translation
// invariance -- the "water-moved" cases that proved important for catching geometry-key / cache bugs.
//
// Migrated off the direct-SCFIterator scaffold (OpenWork D): the facade IS the SALC recipe now.  The
// invariance tolerances need a tightly converged SCF, so each run re-Converges with the tight params the
// scaffold used (the ctor's default convergence is a throwaway warm-up); anchors are unchanged.
#include "gtest/gtest.h"
#include <cmath>

import qchem.Calculation;            // Calculation, CalcOptions, Model, Pol
import qchem.Structure;              // Molecule, Atom
import qchem.SCFIterator;            // SCFParams, EnergyBreakdown
import qchem.Types;                  // Vector3D
import qchem.Math;                   // cos, sin (for the rotation test)
using namespace qchem;

// Tight convergence (== the old M_Sym RunHF/RunDFT SCFParams): drives Δρ/[F,D] to ~1e-9 so the
// symmetric and non-symmetric energies agree to the 1e-6 invariance tolerance.
static const SCFParams tight = {.NMaxIter = 60, .MinΔρ = 1e-7, .MinΔFD = 1e-9, .MinVirial = 1e2,
                                .MinFD = 1e-7, .StartingRelaxRo = 0.5, .MergeTol = 1e-4, .Verbose = false};

static Molecule MakeWater()
{
    // Experimental geometry in BOHR (O-H = 1.809 a0, angle 104.5 deg), C2 axis along z.
    Molecule w;
    w.Insert(new Atom(8, 0, Vector3D<double>(0, 0.0,    0.0)));
    w.Insert(new Atom(1, 0, Vector3D<double>(0, 1.431,  1.107)));
    w.Insert(new Atom(1, 0, Vector3D<double>(0,-1.431,  1.107)));
    return w;
}

// One converged run through the facade.  symmetry=true SALC-blocks the (Cartesian PG) basis; everything
// else is identical between the two, so any energy difference is a SALC-transform error.
static EnergyBreakdown Run(const Molecule& mol, Model model, Pol pol, bool symmetry,
                           Angular angular = Angular::Cartesian)
{
    // The facade's default seed (auto SAD for DFT, CoreGuess for HF) -- the SAD polarized-DFT path is now
    // robust (FittedVxcPol handles the spin-agnostic seed; see M_DFT.WaterPolarizedSAD).  This is an
    // invariance test anyway, so the converged energy is seed-independent.
    Calculation calc(mol, {.basis = "dzvp", .model = model, .pol = pol, .angular = angular, .symmetry = symmetry});
    calc.Converge(tight);                 // re-converge tight for the invariance tolerance
    return calc.EnergyTerms();
}

// Symmetry-adapted Xalpha-DFT water (fitted Coulomb + fitted Vxc) must equal the non-symmetric DFT run:
// the decorator transforms the 3-centre integrals by O, like J/K, building each from the raw basis's
// *cached* 3C (re-entrant integral cache) -- so the raw 3C is shared across irreps.  (xalpha defaults to
// 0.7, == the scaffold's RunDFT alpha; the facade auto-runs DFT with DIIS-from-start.)
static void CheckWaterDFT(Pol pol, double tol)
{
    const Molecule water = MakeWater();
    EnergyBreakdown ebRef = Run(water, Model::Xalpha, pol, false);
    EnergyBreakdown ebSym = Run(water, Model::Xalpha, pol, true);
    EXPECT_NEAR(ebSym.GetTotalEnergy(), ebRef.GetTotalEnergy(), tol) << "DFT symmetric == non-symmetric";
}

TEST(M_Sym, water_DFT_unpolarized) { CheckWaterDFT(Pol::UnPolarized, 1e-5); }
// Polarized (unrestricted) Xalpha: with DIIS driving from the start the non-symmetric and SALC runs
// converge to the SAME stationary point to ~12 digits -- the SALC transform is exact.  (See
// memory project_msym_layout_ub: the old loose 1e-3 was masking a non-convergent run.)
TEST(M_Sym, water_DFT_polarized)   { CheckWaterDFT(Pol::Polarized, 1e-6); }

// Symmetry-adapted water HF (real DZVP basis) must equal the non-symmetric run, both for the closed-shell
// unpolarized and (since water is closed shell) the polarized Hamiltonian -- and the absolute energy /
// virial must be physical.
static void CheckWaterHF(Pol pol)
{
    const Molecule water = MakeWater();
    EnergyBreakdown ebRef = Run(water, Model::HF, pol, false);
    EnergyBreakdown ebSym = Run(water, Model::HF, pol, true);

    EXPECT_NEAR(ebSym.GetTotalEnergy(), ebRef.GetTotalEnergy(), 1e-6) << "symmetric == non-symmetric";
    // physical sanity (water HF/DZVP): total energy near -76 Ha, virial 2+V/K near 0.
    EXPECT_NEAR(ebRef.GetTotalEnergy(), -76.0, 0.5);
    EXPECT_NEAR(ebRef.GetVirial(), -2.0, 0.02) << "virial 2+V/K should be ~0 at the HF minimum";
}

TEST(M_Sym, water_HF_unpolarized) { CheckWaterHF(Pol::UnPolarized); }
TEST(M_Sym, water_HF_polarized)   { CheckWaterHF(Pol::Polarized); }

// SPHERICAL SALC end-to-end (OpenWork A, S3a-S5): the in-house MnD-spherical basis (real solid harmonics;
// water/dzvp carries the O d-shell, so the harmonic path is genuinely exercised) symmetry-adapted must
// equal the un-adapted spherical run.  This drives ExtractAoShells(SphData) + SphericalShellRep through the
// SALC pipeline and the SymmetryAdapt SphData dispatch.  We assert only adapted == un-adapted (the SALC
// invariant); the spherical absolute energy differs from Cartesian by the dropped d s-contaminant (see
// M_DFT.WaterSpherical), so we bound it loosely for physical sanity.
static void CheckWaterHFSpherical(Pol pol)
{
    const Molecule water = MakeWater();
    EnergyBreakdown ebRef = Run(water, Model::HF, pol, false, Angular::Spherical);
    EnergyBreakdown ebSym = Run(water, Model::HF, pol, true,  Angular::Spherical);
    EXPECT_NEAR(ebSym.GetTotalEnergy(), ebRef.GetTotalEnergy(), 1e-6) << "spherical: symmetric == non-symmetric";
    EXPECT_NEAR(ebRef.GetTotalEnergy(), -76.0, 0.5) << "physical sanity (water HF ~ -76 Ha)";
}
TEST(M_Sym, water_HF_spherical_unpolarized) { CheckWaterHFSpherical(Pol::UnPolarized); }
TEST(M_Sym, water_HF_spherical_polarized)   { CheckWaterHFSpherical(Pol::Polarized); }

// Rotate v by Euler (Rz(a) Ry(b) Rx(c)) -- a generic, non-axis-aligned orientation.
static Vector3D<double> Rotate(const Vector3D<double>& v, double a, double b, double c)
{
    double x=v.x, y=v.y, z=v.z;
    double y1=cos(c)*y - sin(c)*z, z1=sin(c)*y + cos(c)*z;                 // Rx(c)
    double x2=cos(b)*x  + sin(b)*z1, z2=-sin(b)*x + cos(b)*z1;             // Ry(b)
    double x3=cos(a)*x2 - sin(a)*y1, y3=sin(a)*x2 + cos(a)*y1;             // Rz(a)
    return Vector3D<double>(x3, y3, z2);
}

// The canonical water, rigidly rotated to a generic orientation and translated off the origin.
static Molecule MakeWaterMoved(double a=0.7, double b=1.1, double c=0.3,
                               Vector3D<double> T=Vector3D<double>(2.0,-3.0,1.5))
{
    Vector3D<double> O(0,0,0), H1(0,1.431,1.107), H2(0,-1.431,1.107);
    Molecule w;
    w.Insert(new Atom(8, 0, Rotate(O ,a,b,c)+T));
    w.Insert(new Atom(1, 0, Rotate(H1,a,b,c)+T));
    w.Insert(new Atom(1, 0, Rotate(H2,a,b,c)+T));
    return w;
}

// Symmetric == non-symmetric HF for a water displaced off the canonical frame, AND equal to the canonical
// -76.0229 (HF is invariant under rigid rotation/translation).  Exercises the whole pipeline -- detection,
// SALCs and the 2-e decorator -- in a generic orientation.  (Relies on the 2-electron integral cache being
// geometry-aware; the canonical run is in the same process.)
static void CheckMovedWaterHF(const Molecule& mol)
{
    EnergyBreakdown ebRef = Run(mol, Model::HF, Pol::UnPolarized, false);
    EnergyBreakdown ebSym = Run(mol, Model::HF, Pol::UnPolarized, true);

    EXPECT_NEAR(ebSym.GetTotalEnergy(), ebRef.GetTotalEnergy(), 1e-6) << "symmetric == non-symmetric";
    EXPECT_NEAR(ebSym.GetTotalEnergy(), -76.022903, 1e-4) << "invariant under rigid rotation/translation";
}

TEST(M_Sym, water_translated)         { CheckMovedWaterHF(MakeWaterMoved(0,0,0,      Vector3D<double>(2.0,-3.0,1.5))); }
TEST(M_Sym, water_rotated)            { CheckMovedWaterHF(MakeWaterMoved(0.7,1.1,0.3, Vector3D<double>(0,0,0)));       }
TEST(M_Sym, water_rotated_translated) { CheckMovedWaterHF(MakeWaterMoved(0.7,1.1,0.3, Vector3D<double>(2.0,-3.0,1.5)));}
