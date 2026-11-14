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
