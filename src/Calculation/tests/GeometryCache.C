// File: Calculation/tests/GeometryCache.C  THE INVARIANT: a geometry iteration builds a new Calculation; no cache outlives its geometry.
//
// doc/IntegralCacheOwnershipPlan.md S2b (lane H, 2026-11-13; the user: the cache must bust between geometry-optimization ITERATIONS, "not
// runs iterations").  The geometry driver steps a surface: GeometryOptimizer's PESCost does R = map.At(q), then surface.ConvergeAt(R), and
// MoleculeSession::ConvergeAt builds a NEW MoleculeCalculation at R (a new basis, a new Hamiltonian, a new Coulomb-fit slot), carrying only
// the previous point's coefficients.  So every integral a point uses, owned by that point or keyed by its exact geometry, is the point's
// own.  If someone later RE-POINTS a live Calculation at a new Structure (or lets any cache outlive its geometry), THIS is the test that
// must catch it: step 2's Coulomb-fit metric, its inverse, the orbital overlap and the orbital x fit 3-centre Coulomb tensor must DIFFER
// from step 1's and EQUAL a cold run's at step 2's geometry.
//
// Driven through MoleculeSession::ConvergeAt, the surface the optimizer steps, not through GeometryOptimizer itself: only a DFT run has
// a Coulomb fit, and LDA's molecular XC route (δ on the Becke mesh) refuses its force today (a MissingGradients() entry), so the
// optimizer's pre-flight would refuse it.  The iteration under test is the surface's step; the optimizer adds only the choice of R.
// "Cold" is a fresh session at step 2's geometry after the first one is gone and the string-keyed global integral cache is EMPTY.
#include "gtest/gtest.h"
#include <cmath>
#include <cstdio>
#include <memory>
#include <optional>
#include <vector>
#include "forward.H"

import qchem.MoleculeCalculation;
import qchem.Hamiltonian.Factory;   // Hamiltonian::Model
import qchem.StructureData;
import qchem.Structure;
import qchem.SCFParams;
import qchem.Types;
import qchem.Blaze;
import qchem.BasisSet;
import qchem.BasisSet.Orbital_DFT_IBS;
import qchem.BasisSet.PairProjection;          // DensePairProjectionIntegrator: one dispatching contraction of the 3-centre record
import qchem.BasisSet.Internal.DB_Cache_RAM;   // the global cache's Clear hooks (tests may cheat): a COLD run in-process
using namespace qchem;

//! What one geometry point holds, read through the point (a friend: the facade hands none of it out).
class GeometryCacheTests
{
public:
    struct Snapshot
    {
        rsmat_t                metric, inverse;   //!< the Coulomb fit's metric and its inverse
        std::vector<rsmat_t>   overlap;           //!< per orbital block
        std::vector<rsmat_t>   repulsion3C;       //!< per orbital block: Σ_f a_f ⟨ij|f⟩ at a fixed a (the dense 3-centre tensor, contracted)
    };
    static Snapshot Take(const MoleculeCalculation& c)
    {
        Snapshot s;
        const auto* slot=c.itsCoulombFit.Read();
        EXPECT_NE(slot, nullptr) << "the point's Coulomb fit slot was never filled";
        if (!slot || !*slot) return s;
        // The Dunlap face (an abstract-to-abstract cross-cast): the metric, its inverse and the charges live there.
        const auto* nonOrtho=dynamic_cast<const BasisSet::FIT_CD_NonOrtho*>(slot->get());
        EXPECT_NE(nonOrtho, nullptr) << "the run's Coulomb fit is not the Dunlap (non-orthogonal) kind";
        if (!nonOrtho) return s;
        const BasisSet::FIT_CD_NonOrtho& fit=*nonOrtho;
        s.metric=fit.Repulsion();
        s.inverse=fit.InvRepulsion();
        rvec_t a(fit.GetNumFunctions());
        for (size_t f=0;f<a.size();++f) a[f]=1.0/(1.0+f);
        for (auto ibs : c.itsBasis->Iterate<BasisSet::Real_OIBS>())
        {
            s.overlap.push_back(ibs->Overlap());
            const auto* dft=dynamic_cast<const BasisSet::Orbital_DFT_IBS<double,double>*>(ibs);   // abstract face to abstract face
            EXPECT_NE(dft, nullptr);
            if (dft) s.repulsion3C.push_back(DensePairProjectionIntegrator<double>(dft->Repulsion3C(fit), fit.Charge()).Adjoint(a));
        }
        return s;
    }
};

namespace
{
void ClearIntegralCache()   //!< as CacheKeyExactness: every operator of the string-keyed global cache
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
template <class M> double MaxAbsDiff(const M& x, const M& y)
{
    EXPECT_EQ(x.rows(), y.rows()); EXPECT_EQ(x.columns(), y.columns());
    double d=0.0;
    for (size_t i=0;i<x.rows();++i) for (size_t j=0;j<x.columns();++j) d=std::max(d, std::fabs(x(i,j)-y(i,j)));
    return d;
}
//! The largest element difference of every matrix in two snapshots, and the largest element (for a relative judgement).
struct Compared {double diff=0.0, scale=0.0;};
Compared Compare(const GeometryCacheTests::Snapshot& x, const GeometryCacheTests::Snapshot& y)
{
    Compared c;
    auto one=[&](const rsmat_t& a, const rsmat_t& b)
    {
        c.diff=std::max(c.diff, MaxAbsDiff(a, b));
        for (size_t i=0;i<a.rows();++i) for (size_t j=0;j<a.columns();++j) c.scale=std::max(c.scale, std::fabs(a(i,j)));
    };
    one(x.metric, y.metric);
    one(x.inverse, y.inverse);
    EXPECT_EQ(x.overlap.size(), y.overlap.size());
    EXPECT_EQ(x.repulsion3C.size(), y.repulsion3C.size());
    for (size_t b=0;b<std::min(x.overlap.size(), y.overlap.size());++b) one(x.overlap[b], y.overlap[b]);
    for (size_t b=0;b<std::min(x.repulsion3C.size(), y.repulsion3C.size());++b) one(x.repulsion3C[b], y.repulsion3C[b]);
    return c;
}
} //anon

TEST(GeometryCache, AGeometryIterationBuildsANewCalculationAndNoCacheOutlivesItsGeometry)
{
    const auto setup=MoleculeSetup::Resolve(StructureData::GetMolecule("H2O"), {.basis="dzvp", .model=Hamiltonian::Model::LDA});
    SCFParams scf;   // the default gate: the claim is about integrals, not the SCF tail
    std::vector<rvec3_t> R1, R2;
    GeometryCacheTests::Snapshot s1, s2, cold;
    ClearIntegralCache();
    {
        MoleculeSession pes(setup, scf);
        R1=pes.Reference().R;
        R2=R1;
        R2.back().y+=0.05;   // one H moves 0.05 bohr: a geometry step
        ASSERT_TRUE(bool(pes.ConvergeAt(std::nullopt, R1)));
        s1=GeometryCacheTests::Take(*pes.Current());
        ASSERT_TRUE(bool(pes.ConvergeAt(std::nullopt, R2)));   // the SECOND iteration of one surface
        s2=GeometryCacheTests::Take(*pes.Current());
    }   // the surface and its points are gone before the cache is cleared (its Get hands out references)
    ClearIntegralCache();
    {
        MoleculeSession fresh(setup, scf);
        ASSERT_TRUE(bool(fresh.ConvergeAt(std::nullopt, R2)));
        cold=GeometryCacheTests::Take(*fresh.Current());
    }
    const Compared moved=Compare(s2, s1), stale=Compare(s2, cold);
    std::printf("  [geometry cache] step 2 vs step 1: max|Δ| %.3e;  step 2 vs a cold run at step 2: max|Δ| %.3e  (largest element %.3e)\n",
                moved.diff, stale.diff, stale.scale);
    EXPECT_GT(moved.diff, 1e-6) << "step 2's fit metric, overlap and 3-centre tensor must move with the atom: a cache outlived its geometry";
    EXPECT_LE(stale.diff, 1e-13*stale.scale) << "step 2 must equal a cold run at its geometry: something carried step 1's integrals";
}
