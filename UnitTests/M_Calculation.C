// File: UnitTests/M_Calculation.C  Tests for the qchem::Calculation front-door facade.
//
// The facade is meant to BE the production recipe -- so its job here is to reproduce, through the
// public front door alone, the converged HF/dzvp water energy that M_HF_U checks via the QchemTester
// scaffold.  Same geometry, same basis, same anchor.  Also a quick smoke of the ScalarFunction
// sampling surface (Density/HOMO/Orbital), the win the API-ergonomics review highlighted.
#include "gtest/gtest.h"
#include <cmath>

import qchem.Calculation;
import qchem.Structure;
import qchem.Types;        // Vector3D

using qchem::Calculation;

static Molecule MakeWater()      // experimental geometry in BOHR, C2 axis along z (== M_HF_U)
{
    Molecule w;
    w.Insert(new Atom(8, 0, Vector3D<double>(0,  0.0,   0.0)));
    w.Insert(new Atom(1, 0, Vector3D<double>(0,  1.431, 1.107)));
    w.Insert(new Atom(1, 0, Vector3D<double>(0, -1.431, 1.107)));
    return w;
}

// Regression anchor: identical converged HF/dzvp total energy to M_HF_U_Water (-76.022903).
TEST(M_Calculation, WaterEnergy)
{
    Molecule water = MakeWater();
    Calculation calc(water, {.basis = "dzvp"});            // build + converge in one line

    const double E_ref = -76.022903;
    EXPECT_LT(std::fabs((E_ref - calc.Energy()) / E_ref), 1e-5);
    EXPECT_GT(calc.IterationCount(), 0u);
}

// {.symmetry=true} SALC-blocks the (Cartesian PG) basis + does global aufbau across irreps.  The
// converged total energy must match the un-blocked run -- the M_Sym invariant, now through the facade.
TEST(M_Calculation, WaterSymmetry)
{
    Calculation calc(MakeWater(), {.basis = "dzvp", .symmetry = true});

    const double E_ref = -76.022903;                       // identical to the un-blocked WaterEnergy anchor
    EXPECT_LT(std::fabs((E_ref - calc.Energy()) / E_ref), 1e-5);
    EXPECT_GT(calc.IterationCount(), 0u);
}

// The caller's Molecule is deep-copied: it survives being passed in and may be used afterwards.
TEST(M_Calculation, OwnsItsOwnStructure)
{
    Calculation calc(MakeWater(), {.basis = "dzvp"});      // temporary destroyed at the semicolon
    EXPECT_EQ(calc.GetStructure().GetNumAtoms(), 3u);
    EXPECT_NEAR(calc.GetStructure().GetNumElectrons(), 10.0, 1e-12);
}

// Density and every occupied MO are sampleable ScalarFunction<double>s through one interface.
TEST(M_Calculation, SamplingSurface)
{
    Calculation calc(MakeWater(), {.basis = "dzvp"});

    // Water has 10 electrons -> 5 doubly-occupied MOs in the unpolarized wave function.
    EXPECT_EQ(calc.NumOccupied(), 5u);

    const Vector3D<double> rO(0, 0, 0);                    // at the oxygen nucleus
    EXPECT_GT(calc.Density()(rO), 0.0);                    // rho > 0 where the atoms are
    EXPECT_GT(std::fabs(calc.HOMO()(rO)) + std::fabs(calc.Orbital(0)(rO)), 0.0);
}
