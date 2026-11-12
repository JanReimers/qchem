// File: Calculation/tests/CacheKeyExactness.C  Two geometries 3e-6 bohr apart, in either order in ONE process, each give their own
// fresh-process energy (lane H, 2026-11-12; OpenWork §3 "a cache key is exact", found by lane A).
//
// The process-wide integral cache keyed every basis and structure on a STREAMED double (6 significant digits), so the second of two
// H₂O₂ geometries 3e-6 bohr apart was served the first one's integrals: 2.8e-5 Ha wrong, and a BFGS relaxation that revisited within
// 1e-6 bohr took a spurious minimum.  "Fresh" here is a run with the string-keyed integral cache EMPTY (the earlier facades are gone
// before it is cleared: its Get hands out references).  The primitive caches (Cache2/3) key on exact doubles and are left alone.
#include "gtest/gtest.h"
#include <cmath>
#include <cstdio>

import qchem.MoleculeCalculation;
import qchem.Hamiltonian.Factory;   // Hamiltonian::Model
import qchem.Structure;
import qchem.SCFParams;
import qchem.Types;
import qchem.BasisSet;
import qchem.BasisSet.Internal.DB_Cache_RAM;   // the concrete cache's Clear hooks (tests may cheat): a FRESH evaluation in-process
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
    MoleculeCalculation calc(MoleculeSetup::Resolve(m, {.basis="dzvp", .model=Hamiltonian::Model::HF}), p);
    EXPECT_TRUE(calc.IsConverged());
    return calc.Energy();
}
//! Empty the string-keyed integral cache: every operator of every map.
void ClearIntegralCache()
{
    using C=BasisSet::IntegralsCache_Base;
    auto& ram=dynamic_cast<BasisSet::IntegralsCache_RAM<double>&>(BasisSet::theCache<double>());
    for (auto op : {C::I1C::Charge, C::I1C::Normalization}) ram.Clear(op);
    for (auto op : {C::I2C::Overlap, C::I2C::Repulsion, C::I2C::Kinetic, C::I2C::Nuclear, C::I2C::RestMass, C::I2C::InvOverlap,
                    C::I2C::InvRepulsion}) ram.Clear(op);
    for (auto op : {C::I2n::Nuclear, C::I2n::LocalPP, C::I2n::LocalPPLong, C::I2n::LocalPPShort, C::I2n::SeparablePP}) ram.Clear(op);
    for (auto op : {C::I2x::Kinetic, C::I2x::Repulsion, C::I2x::Overlap}) ram.Clear(op);
    for (auto op : {C::I3C::Overlap, C::I3C::Repulsion}) ram.Clear(op);
    for (auto op : {C::I4C::Direct, C::I4C::Exchange}) ram.Clear(op);
}
} //anon

TEST(CacheKeyExactness, TwoGeometries3e6BohrApartEachGiveTheirFreshEnergyInEitherOrder)
{
    const double dx=3e-6;
    ClearIntegralCache();
    const double eA=H2O2(0.0);         // fresh
    const double eBafterA=H2O2(dx);    // the cache holds A's integrals
    ClearIntegralCache();
    const double eB=H2O2(dx);          // fresh
    const double eAafterB=H2O2(0.0);   // the cache holds B's integrals
    std::printf("  [cache key] E(A) %.12f  E(B) %.12f  E(B|A)-E(B) %.2e  E(A|B)-E(A) %.2e\n", eA, eB, eBafterA-eB, eAafterB-eA);
    EXPECT_NEAR(eBafterA, eB, 1e-10) << "B after A was served A's integrals (the key printed 6 significant digits)";
    EXPECT_NEAR(eAafterB, eA, 1e-10) << "A after B was served B's integrals";
    EXPECT_GT(std::fabs(eA-eB), 1e-12) << "the two geometries must be different points for the gate to mean anything";
}
