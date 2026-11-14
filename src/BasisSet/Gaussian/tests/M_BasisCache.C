// File: M_BasisCache.C  The parent-made basis cache and the partner's instance identity (doc/IntegralCacheOwnershipPlan.md S3a, §3b).
//
// The orbital parent's cache keys ⟨ab|c⟩ by (the block, the fit's INSTANCE serial), never by the fit's address or content.  The §3b gate:
//   1. build a fit basis, fill its 3-centre table on an orbital block, destroy the fit;
//   2. the orbital cache frees EXACTLY those entries;
//   3. build a second fit at the SAME address (placement-new into the freed storage forces it);
//   4. it gets a different serial, MISSES, and fills fresh.
// The address-reuse half is the bug class the token exists for: in June 2026 a density keyed by its pointer was served a freed one's
// matrices at the same address (M_DFT_Water diverged at iteration 5, allocator-dependent).
#include "gtest/gtest.h"
#include <cstdint>
#include <memory>
#include <new>
#include <string>
#include <cmath>
#include <blaze/math/expressions/DMatDMatEqualExpr.h>   // op== inside gtest (as DBCache.C)

import qchem.BasisSet.Gaussian.Point.Factory;                 // Factory, BasisSetData, Engine, Angular
import qchem.BasisSet.Gaussian.PG_Cart;                       // EFit_IBS: the concrete Coulomb fit (placement-new needs the type)
import qchem.BasisSet.Gaussian.Point.Readers.Gaussian94;
import qchem.BasisSet.Gaussian.Point.BasisFiles;
import qchem.BasisSet.Internal.BasisSetImp;                   // the parent's BlockCache() (tests may cheat)
import qchem.BasisSet;
import qchem.BasisSet.Orbital_DFT_IBS;
import qchem.Structure;
import qchem.Types;
import qchem.Mesh;                                            // MeshParams (the δ fit's XC mesh)
import qchem.Symmetry;                                        // SymmetryOrder: the old δ table key
using namespace qchem;
using namespace qchem::BasisSet::Gaussian;

TEST(M_BasisCache, ATokenAtAReusedAddressIsANewInstance)
{
    struct Holder {InstanceToken token;};
    alignas(Holder) unsigned char storage[sizeof(Holder)];
    Holder* first=new (storage) Holder;
    const uint64_t s1=first->token.Serial();
    const auto alive1=first->token.Alive();
    EXPECT_FALSE(alive1.expired());
    first->~Holder();
    EXPECT_TRUE(alive1.expired()) << "the weak handle must see its owner die";
    Holder* second=new (storage) Holder;
    EXPECT_EQ(static_cast<void*>(first), static_cast<void*>(second)) << "placement-new did not reuse the address: the gate would be vacuous";
    EXPECT_NE(second->token.Serial(), s1) << "a new object at a reused address took the dead one's serial";
    const Holder copy(*second);
    EXPECT_NE(copy.token.Serial(), second->token.Serial()) << "a copy is a new instance";
    second->~Holder();
}

TEST(M_BasisCache, ADeadFitsTablesAreFreedAndANewFitAtItsAddressMisses)
{
    Molecule mol;
    mol.Insert(Site{8, rvec3_t( 0.0,     0.0,      0.0)});
    mol.Insert(Site{1, rvec3_t( 1.43042, 1.10715,  0.0)});
    mol.Insert(Site{1, rvec3_t(-1.43042, 1.10715,  0.0)});
    std::unique_ptr<BasisSet::Real_BS> bs(Factory(BasisSetData::DZVP, &mol, Engine::MnD, Angular::Cartesian));
    const auto& parent=dynamic_cast<const BasisSet::BasisSetImp<double>&>(*bs);   // a cheat: the cache the set MADE
    const BasisSet::BasisCache& cache=*parent.BlockCache();
    const BasisSet::Real_DFT_OIBS* ibs=nullptr;
    for (auto b : bs->Iterate<BasisSet::Real_DFT_OIBS>()) {ibs=b; break;}
    ASSERT_NE(ibs, nullptr);

    using Fit=PG_Cart::EFit_IBS;
    auto makeFit=[&](void* where)
    {
        Gaussian94Reader reader(BasisFile("A1_coul.bsd"));
        return new (where) Fit(std::make_shared<BasisSet::BasisCache>(), &reader, &mol, "A1_coul");
    };
    alignas(Fit) unsigned char storage[sizeof(Fit)];

    const size_t n0=cache.NumThreeCentre();
    Fit* fit1=makeFit(storage);
    const uint64_t serial1=fit1->Instance().Serial();
    const PairProjection<double>& R1=ibs->Repulsion3C(*fit1);
    const auto dense1=R1.dense;                                        // a copy: the stored one dies with its partner
    EXPECT_EQ(&R1, &ibs->Repulsion3C(*fit1)) << "the same fit asked twice must be a HIT";
    EXPECT_EQ(cache.NumThreeCentre(), n0+1);

    fit1->~Fit();
    EXPECT_EQ(cache.EvictDeadPartners(), 1u) << "exactly the dead fit's table is freed";
    EXPECT_EQ(cache.NumThreeCentre(), n0);
    EXPECT_EQ(cache.EvictDeadPartners(), 0u) << "a second sweep finds nothing";

    Fit* fit2=makeFit(storage);
    EXPECT_EQ(static_cast<void*>(fit1), static_cast<void*>(fit2)) << "placement-new did not reuse the address: the gate would be vacuous";
    EXPECT_NE(fit2->Instance().Serial(), serial1);
    const PairProjection<double>& R2=ibs->Repulsion3C(*fit2);           // a MISS: filled fresh
    EXPECT_EQ(cache.NumThreeCentre(), n0+1);
    EXPECT_EQ(R2.dense, dense1) << "the same fit content rebuilt gives the same numbers";
    fit2->~Fit();
}

// S3c, the latent defect the δ basis's own tables carried: they were keyed by the block's SPATIAL symmetry, so two blocks of equal
// symmetry (two basis sets, or one rebuilt at a new geometry, against a surviving fit) got the FIRST block's Φ.  Now the fit keys
// them by the block's InstanceToken: block 2 gets its own table, the one a fresh fit builds for it.
TEST(M_BasisCache, TwoBlocksOfEqualSymmetryAskedOfOneDeltaFitGetTheirOwnTables)
{
    auto water=[](double dy)
    {
        Molecule m;
        m.Insert(Site{8, rvec3_t( 0.0,     0.0,      0.0)});
        m.Insert(Site{1, rvec3_t( 1.43042, 1.10715,  0.0)});
        m.Insert(Site{1, rvec3_t(-1.43042, 1.10715+dy, 0.0)});
        return m;
    };
    const Molecule mol1=water(0.0), mol2=water(0.1);   // one H moved 0.1 bohr: the same symmetry, different functions
    std::unique_ptr<BasisSet::Real_BS> bs1(Factory(BasisSetData::DZVP, &mol1, Engine::MnD, Angular::Cartesian));
    std::unique_ptr<BasisSet::Real_BS> bs2(Factory(BasisSetData::DZVP, &mol2, Engine::MnD, Angular::Cartesian));
    auto firstBlock=[](const BasisSet::Real_BS& bs)
    {
        const BasisSet::Real_DFT_OIBS* b=nullptr;
        for (auto x : bs.Iterate<BasisSet::Real_DFT_OIBS>()) {b=x; break;}
        return b;
    };
    const auto* b1=firstBlock(*bs1);
    const auto* b2=firstBlock(*bs2);
    ASSERT_TRUE(b1 && b2);
    const SymmetryOrder less;
    ASSERT_FALSE(less(b1->GetSymt(), b2->GetSymt()) || less(b2->GetSymt(), b1->GetSymt()))
        << "the two blocks must share a spatial symmetry, or the gate is vacuous";

    const qcMesh::MeshParams mesh{.radial=qcMesh::RadialKind::MHL, .nRadial=40, .mhl_m=3, .mhl_alpha=2.0,
                                  .angular=qcMesh::AngularKind::Lebedev, .angularDegree=17, .beckeOrder=2};
    BasisSet::FitQuadrature quad=BasisSet::MakeXCQuadrature(mol1, mesh, nullptr);
    std::unique_ptr<const BasisSet::rFIT_SF_ABS> delta(bs1->CreateVxcFitBasisSet(&mol1, mesh, BasisSet::VxcFit::Delta, &quad));
    std::unique_ptr<const BasisSet::rFIT_SF_ABS> fresh(bs1->CreateVxcFitBasisSet(&mol1, mesh, BasisSet::VxcFit::Delta, &quad));

    const PairProjection<double>& r1=delta->Overlap3C(*b1);
    const PairProjection<double>& r2=delta->Overlap3C(*b2);
    const PairProjection<double>& f2=fresh->Overlap3C(*b2);
    EXPECT_NE(&r1, &r2) << "block 2 was served block 1's record";
    // Compare the records by what they DO: the adjoint Φ† diag(w c) Φ for one coefficient vector over the points.
    const size_t np=delta->GetNumFunctions();
    rvec_t c(np);
    for (size_t g=0; g<np; g++) c[g]=1.0+0.37*std::sin(0.11*double(g));
    const auto A1=r1.applyRawAdjoint(c), A2=r2.applyRawAdjoint(c), F2=f2.applyRawAdjoint(c);
    ASSERT_EQ(A2.rows(), F2.rows()); ASSERT_EQ(A2.rows(), A1.rows());
    double moved=0.0, stale=0.0;
    for (size_t i=0; i<A2.rows(); i++)
        for (size_t j=0; j<A2.columns(); j++)
        {
            moved=std::max(moved, std::fabs(A2(i,j)-A1(i,j)));
            stale=std::max(stale, std::fabs(A2(i,j)-F2(i,j)));
        }
    EXPECT_GT(moved, 1e-3) << "block 2's record must differ from block 1's (its H moved)";
    EXPECT_EQ(stale, 0.0) << "block 2's record must be the one a fresh fit builds for it";
    EXPECT_EQ(&r1, &delta->Overlap3C(*b1)) << "asked again, block 1 is a HIT";
}
