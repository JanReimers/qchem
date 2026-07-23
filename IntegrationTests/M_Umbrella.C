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
    water.Insert(new Atom(8, 0, Vector3D<double>(0,  0.0,   0.0)));
    water.Insert(new Atom(1, 0, Vector3D<double>(0,  1.431, 1.107)));
    water.Insert(new Atom(1, 0, Vector3D<double>(0, -1.431, 1.107)));

    // qchem::Calculation + CalcOptions designated init  (qchem.Calculation, re-exported)
    qchem::Calculation calc(water, {.basis = "dzvp"});
    calc.Converge(SCFParams{.NMaxIter = 20});            // SCFParams visible (via qchem.SCFIterator)

    const double E_ref = -76.022903;
    EXPECT_LT(std::fabs((E_ref - calc.Energy()) / E_ref), 1e-5);

    // ScalarFunction<double> sampling  (qchem.ScalarFunction, re-exported)
    const ScalarFunction<double>& rho = calc.Density();
    EXPECT_GT(rho(Vector3D<double>(0, 0, 0)), 0.0);
}
