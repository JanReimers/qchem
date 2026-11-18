// File: Calculation/tests/CacheKeyExactness.C  Two geometries 3e-6 bohr apart, in either order in ONE process, each give their own
// fresh-process energy (lane H, 2026-11-12; OpenWork §3 "a cache key is exact", found by lane A).
//
// The process-wide integral cache keyed every basis and structure on a STREAMED double (6 significant digits), so the second of two
// H₂O₂ geometries 3e-6 bohr apart was served the first one's integrals: 2.8e-5 Ha wrong, and a BFGS relaxation that revisited within
// 1e-6 bohr took a spurious minimum.  Since S6 there is no process-wide cache: each run owns its caches (doc/IntegralCacheOwnershipPlan.md),
// so every run is fresh, and the gate guards that no future cache brings the cross-geometry serving back.
#include "gtest/gtest.h"
#include <cmath>
#include <cstdio>

import qchem.MoleculeCalculation;
import qchem.Hamiltonian.Factory;   // Hamiltonian::Model
import qchem.Structure;
import qchem.SCFParams;
import qchem.Types;
import qchem.BasisSet;
using namespace qchem;

namespace
{
//! H₂O₂ HF/DZVP with the first O at x = -1.322808 + \a dx, converged tight; its energy.  The facade dies before returning.
double H2O2(double dx)
{
    Molecule m;
    m.Insert(Site{8, rvec3_t(-1.322808+dx, 0.0, 0.0)});
    m.Insert(Site{8, rvec3_t( 1.322808, 0.0, 0.0)});
    m.Insert(Site{1, rvec3_t(-0.949557, 1.520749, 0.878005)});
    m.Insert(Site{1, rvec3_t( 0.949557,-1.520749, 0.878005)});
    SCFParams p; p.NMaxIter=100; p.MinΔρPerElectron=1e-9;
    MoleculeCalculation calc(MoleculeSetup::Resolve(m, {.model=Hamiltonian::Model::HF, .molecule={.basis="dzvp"}}), p);
    EXPECT_TRUE(calc.IsConverged());
    return calc.Energy();
}
//! (S6: every run builds its own basis set and so its own caches; there is no process-wide cache left to clear, so every run here is
//! fresh by construction.  The gate stays: a second geometry must never be served the first's integrals, whatever caches them.)
} //anon

TEST(CacheKeyExactness, TwoGeometries3e6BohrApartEachGiveTheirFreshEnergyInEitherOrder)
{
    const double dx=3e-6;
    const double eA=H2O2(0.0);         // fresh
    const double eBafterA=H2O2(dx);    // the cache holds A's integrals
    const double eB=H2O2(dx);          // fresh
    const double eAafterB=H2O2(0.0);   // the cache holds B's integrals
    std::printf("  [cache key] E(A) %.12f  E(B) %.12f  E(B|A)-E(B) %.2e  E(A|B)-E(A) %.2e\n", eA, eB, eBafterA-eB, eAafterB-eA);
    EXPECT_NEAR(eBafterA, eB, 1e-10) << "B after A was served A's integrals (the key printed 6 significant digits)";
    EXPECT_NEAR(eAafterB, eA, 1e-10) << "A after B was served B's integrals";
    EXPECT_GT(std::fabs(eA-eB), 1e-12) << "the two geometries must be different points for the gate to mean anything";
}
