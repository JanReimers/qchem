// File: UnitTests/M_Umbrella.C  Proves `import qchem;` alone exposes the molecule front door.
//
// Intentionally imports ONLY the umbrella -- no granular qchem.* imports.  If any type below stops
// being re-exported, this TU fails to compile, which is the whole point: the umbrella is the single
// consumer surface and this pins it.
#include "gtest/gtest.h"
#include <cmath>

import qchem;   // <-- the umbrella, nothing else
using namespace qchem;

TEST(M_Umbrella, FrontDoorViaUmbrella)
{
    // Molecule / Atom / Vector3D  (qchem.Structure + qchem.Types, re-exported)
    Molecule water;
    water.Insert(Site{8, Vector3D<double>(0,  0.0,   0.0)});
    water.Insert(Site{1, Vector3D<double>(0,  1.431, 1.107)});
    water.Insert(Site{1, Vector3D<double>(0, -1.431, 1.107)});

    // qchem::MoleculeCalculation + RunOptions designated init  (qchem.MoleculeCalculation, re-exported)
    auto calc=qchem::MoleculeCalculation::Prepare(MoleculeSetup::Resolve(water, {.molecule={.basis="dzvp"}})).TakeValue();
    ASSERT_TRUE(calc->Converge([](SCFParams& p){ p.NMaxIter = 20; }));   // SCFParams visible; a one-field edit (R10); the Outcome door

    const double E_ref = -76.022903;
    EXPECT_LT(std::fabs((E_ref - calc->Energy()) / E_ref), 1e-5);

    // ScalarFunction<double> sampling  (qchem.ScalarFunction, re-exported)
    const ScalarFunction<double>& rho = calc->Density();
    EXPECT_GT(rho(Vector3D<double>(0, 0, 0)), 0.0);
}
