// File: UnitTests/A_Cache4.C  Unit test the Atom Cach4 system for storing 4 index radial Slater integrals.
#include "gtest/gtest.h"
#include <gmock/gmock.h>
#include <iostream>
#include <map>
#include <string>
#include <memory>
#include "forward.H"
#include <nlohmann/json.hpp>
using std::cout;
using std::endl;

import qchem.BasisSet.Radial.Evaluators; 
import qchem.BasisSet.Radial.Evaluators.Gaussian.IBS; 
import qchem.BasisSet.Radial.Evaluators.Internal.Grouper;
import qchem.BasisSet.Radial.Evaluators.Internal.ExponentialEvaluator;
import qchem.BasisSet.Radial.Factory;
import qchem.BasisSet.Internal.Orbital_ERI4_IBS;  // tests may cheat: the ERI4 substrate face
import qchem.BasisSet.Internal.BasisSetImp;     // the parent's BlockCache() (tests may cheat): its Rk tables
import qchem.BasisSet.Internal.Cache4;         // Cache4::Lookups()/Inserts()
using namespace qchem;

using namespace qchem::BasisSet::Radial;
using namespace qchem::BasisSet::Radial::Evaluators;
using enum BasisSetAccuracy;
using GCache4=BasisSet::Radial::Evaluators::Gaussian::Gaussian_Cache4;

// The tiny fixed test basis the removed BasisSetAccuracy::N3 used to build (3 even-tempered exponents {0.5,1,2}
// / BSpline N=3).  N3/N5 were removed from the production pool (ill-conditioned); this cache test wants a small
// deterministic basis, so build the equivalent JSON directly (test-local).
static BasisSet::Real_BS* N3Basis(Type t, size_t Z)
{
    nlohmann::json js = (t==Type::BSpline6||t==Type::BSpliner6)
        ? nlohmann::json{{"type",t},{"rmin",0.5},{"rmax",2},{"N",3}}
        : nlohmann::json{{"type",t},{"ltrim",0},{"exponents",{0.5,1,2.0}}};
    return Factory(js, Z);
}

class Cache4Tests : public ::testing::Test
{
public:
    Cache4Tests() : bs1(0), bs2(0)
    {

    }

    const ExponentGrouper& GetGrouper(const GCache4* gc)
    {
        assert(gc);
        return gc->grouper; //friend access.
    }

    const std::vector<size_t>& maxls(const GCache4* gc)
    {
        return GetGrouper(gc).maxls;
    }

    const std::vector<size_t>& es_indices(const ExponentialEvaluator* e)
    {
        return e->es_indices;
    }

    void Init(size_t Z, Type type1, Type type2)
    {
        delete bs1;
        delete bs2;
        bs1=N3Basis(type1,Z);
        bs2=N3Basis(type2,Z);
    }
    void InitBSpline6(size_t Z) {Init(Z,Type::BSpline6,Type::BSpline6);}
    void InitGaussian(size_t Z) {Init(Z,Type::Gaussian,Type::Gaussian);}
    void InitSlater  (size_t Z) {Init(Z,Type::Slater  ,Type::Slater);}

    // The ERI4 store is canonical-only (doc/ERI4Rework.md §5.2): only the a<=b block (by the block's index in its parent's cache = the
    // parent's iteration order) is built/stored, and the partner request throws.  So fetch only the canonical block, check the second
    // basis computes the same numbers (in its OWN cache since S3a: no cross-basis object sharing), verify J(a,b)=J(b,a)^T against the
    // UNCACHED MakeDirect of the partner, and assert the non-canonical request is forbidden.
    void TestDirect(double eps, bool testTranspose=true)
    {
        using BasisSet::Real_ERI4_OIBS;
        auto ibs21=bs2->Iterate<Real_ERI4_OIBS>().begin();
        size_t i=0;
        for (auto ibs11:bs1->Iterate<Real_ERI4_OIBS>())
        {
            auto ibs22=bs2->Iterate<Real_ERI4_OIBS>().begin();
            size_t j=0;
            for (auto ibs12:bs1->Iterate<Real_ERI4_OIBS>())
            {
                if (i<=j)
                {
                    const ERI4& J1=ibs11->Direct(*ibs12);
                    const ERI4& J2=(*ibs21)->Direct(**ibs22);
                    EXPECT_EQ(J1,J2);
                    if (testTranspose)
                        EXPECT_NEAR(fnorm(J1,ibs12->MakeDirect(*ibs11).Transpose()),0.0,eps);
                }
                else
                    EXPECT_ANY_THROW(ibs11->Direct(*ibs12));   // non-canonical request forbidden
                ++ibs22; ++j;
            }
            ++ibs21; ++i;
        }

    }
    void TestExchange(double eps, bool testTranspose=true)
    {
        using BasisSet::Real_ERI4_OIBS;
        auto ibs21=bs2->Iterate<Real_ERI4_OIBS>().begin();
        size_t i=0;
        for (auto ibs11:bs1->Iterate<Real_ERI4_OIBS>())
        {
            auto ibs22=bs2->Iterate<Real_ERI4_OIBS>().begin();
            size_t j=0;
            for (auto ibs12:bs1->Iterate<Real_ERI4_OIBS>())
            {
                if (i<=j)
                {
                    const ERI4& K1=ibs11->Exchange(*ibs12);
                    const ERI4& K2=(*ibs21)->Exchange(**ibs22);
                    EXPECT_EQ(K1,K2);
                    if (testTranspose)
                        EXPECT_NEAR(fnorm(K1,ibs12->MakeExchange(*ibs11).Transpose()),0.0,eps);
                }
                else
                    EXPECT_ANY_THROW(ibs11->Exchange(*ibs12));   // non-canonical request forbidden
                ++ibs22; ++j;
            }
            ++ibs21; ++i;
        }

    }
    
    void TestDirect(size_t Z,nlohmann::json js)
    {
        bs1=Factory(js,Z);
        using BasisSet::Real_ERI4_OIBS;
        for (auto ibs:bs1->Iterate<Real_ERI4_OIBS>())
        {
            // std::cout << *ibs;
            ibs->Direct(*ibs);
        }
        delete bs1;
    }
 
    

    // Trigger every 4-index Rk Create() for a basis by evaluating Direct+Exchange over all shell pairs (into the basis's own Rk table).
    void ExerciseAllRkOf(const BasisSet::Real_BS& bs)
    {
        using BasisSet::Real_ERI4_OIBS;
        // Only the canonical (a<=b) block of each pair may be built (the ERI4 cache is canonical-only); the
        // partner exercises the SAME underlying Rk integrals, so canonical pairs still fully warm the pool.
        size_t ii=0;
        for (auto i:bs.Iterate<Real_ERI4_OIBS>())
        {
            size_t jj=0;
            for (auto j:bs.Iterate<Real_ERI4_OIBS>())
            {
                if (ii<=jj)   // canonical: the parent's iteration order IS the cache-index order
                {
                    i->Direct  (*j);
                    i->Exchange(*j);
                }
                jj++;
            }
            ii++;
        }
    }

    // S5a: an atom basis OWNS its Rk table (its parent's BasisCache, doc/IntegralCacheOwnershipPlan.md §1).  Its l blocks share it, so
    // the within-atom reuse stays (the Rk of one exponent quartet serves every l pair), and a second atom gets its own table: the
    // cross-element warm-up the process-wide pool gave is dropped on purpose (user 2026-11-12: "not a very important capability").
    void TestOneTablePerAtom(Type type, const std::string& radialType, size_t Za, size_t Zb)
    {
        std::unique_ptr<BasisSet::Real_BS> a(Factory(Low,type,Za)), b(Factory(Low,type,Zb));
        auto tableOf=[&](const BasisSet::Real_BS& bs) -> const Cache4&
            {return dynamic_cast<const BasisSet::BasisSetImp<double>&>(bs).BlockCache()->GetCache4(radialType);};   // tests may cheat
        const Cache4& ca=tableOf(*a);
        const Cache4& cb=tableOf(*b);
        EXPECT_NE(&ca, &cb) << "two atom basis sets share one Rk table";
        ExerciseAllRkOf(*a);
        const size_t inserts=ca.Inserts(), lookups=ca.Lookups();
        const double reuse=lookups ? 100.0*(1.0-double(inserts)/double(lookups)) : 0.0;
        std::cout << "[S5a] " << radialType << " Z" << Za << ": lookups=" << lookups << " inserts=" << inserts << " reuse=" << reuse << "%" << std::endl;
        EXPECT_GT(lookups, 10000u);                                         // the atom really did substantial Rk work...
        EXPECT_GT(reuse, 50.0) << "the l blocks of one atom must share the Rk of a common exponent quartet";
        EXPECT_EQ(cb.Lookups(), 0u) << "the other atom's table saw none of this atom's work";
    }

    BasisSet::Real_BS *bs1,*bs2;
};

// Do this test first in order to force RkEngine to upgrade itself as l increases.
TEST_F(Cache4Tests,HF2_SG_Reentry)
{
    delete bs1;
    nlohmann::json js={{"N", 5}, {"emin", 0.25}, {"emax", 4}};
    js["type"]=Type::Gaussian;
    for (size_t Z:{2,18,36,86})
        TestDirect(Z,js);

        // Now shift the exponents around
    js["emin"]=0.2;
    js["emax"]=3.0;
    for (size_t Z:{2,18,36,86})
        TestDirect(Z,js);
}
TEST_F(Cache4Tests,HF2_SL_Reentry)
{
    delete bs1;
    nlohmann::json js={{"N", 5}, {"emin", 0.25}, {"emax", 4}};
    js["type"]=Type::Slater;
    for (size_t Z:{2,18,36,86})
        TestDirect(Z,js);

        // Now shift the exponents around
    js["emin"]=0.2;
    js["emax"]=3.0;
    for (size_t Z:{2,18,36,86})
        TestDirect(Z,js);
}

TEST_F(Cache4Tests,HF2_BS_Reentry)
{
    delete bs1;
    nlohmann::json js={{"N", 5}, {"rmin", 0.25}, {"rmax", 4}};
    js["type"]=Type::BSpline6;
    for (size_t Z:{2,18,36,86})
        TestDirect(Z,js);

        // Now shift the exponents around
    js["rmin"]=0.2;
    js["rmax"]=3.0;
    for (size_t Z:{2,18,36,86})
        TestDirect(Z,js);
}

// Re-entry with DIFFERENT exponents must stay distinct.  SG/SL share ONE coarse "SG"/"SL" Cache4
// (RadialType()==Name(), carries no exponents), so correctness rests entirely on the ExponentGrouper
// keying each distinct exponent VALUE to its own index (Slater/Gaussian Evaluator.C use
// grouper.unique_esv[i] directly in the integral).  If that keying ever regressed to a per-basis index,
// the second basis would alias onto the first's cached table and fnorm(JA,JB) would collapse to 0.
// These guard that the coarse RadialType key is SAFE for exponentials.  BSpline gets the same treatment
// now (RadialType()==Name(), one "BSpline<6>" Cache4 for all grids) -- see HF2_BS_GridKeyed, which relies
// on the lossless SplineGrouper key instead of an exponent-value key.
TEST_F(Cache4Tests,HF2_SG_ExponentKeyed)
{
    delete bs1;
    nlohmann::json jsA={{"N", 5}, {"emin", 0.25}, {"emax", 4.0}}; jsA["type"]=Type::Gaussian;
    nlohmann::json jsB={{"N", 5}, {"emin", 0.20}, {"emax", 3.0}}; jsB["type"]=Type::Gaussian;
    BasisSet::Real_BS* bsA=Factory(jsA,2);
    BasisSet::Real_BS* bsB=Factory(jsB,2);   // registers into the SAME "SG" Cache4 as bsA
    using BasisSet::Real_ERI4_OIBS;
    auto a=bsA->Iterate<Real_ERI4_OIBS>().begin();
    auto b=bsB->Iterate<Real_ERI4_OIBS>().begin();
    const ERI4& JA=(*a)->Direct(**a);
    const ERI4& JB=(*b)->Direct(**b);
    EXPECT_GT(fnorm(JA,JB),1e-3);            // distinct exponents -> distinct integrals (no aliasing)
    delete bsA;
    delete bsB;
}
TEST_F(Cache4Tests,HF2_SL_ExponentKeyed)
{
    delete bs1;
    nlohmann::json jsA={{"N", 5}, {"emin", 0.25}, {"emax", 4.0}}; jsA["type"]=Type::Slater;
    nlohmann::json jsB={{"N", 5}, {"emin", 0.20}, {"emax", 3.0}}; jsB["type"]=Type::Slater;
    BasisSet::Real_BS* bsA=Factory(jsA,2);
    BasisSet::Real_BS* bsB=Factory(jsB,2);   // registers into the SAME "SL" Cache4 as bsA
    using BasisSet::Real_ERI4_OIBS;
    auto a=bsA->Iterate<Real_ERI4_OIBS>().begin();
    auto b=bsB->Iterate<Real_ERI4_OIBS>().begin();
    const ERI4& JA=(*a)->Direct(**a);
    const ERI4& JB=(*b)->Direct(**b);
    EXPECT_GT(fnorm(JA,JB),1e-3);            // distinct exponents -> distinct integrals (no aliasing)
    delete bsA;
    delete bsB;
}

// The BSpline analogue of the SG/SL exponent-keyed guard.  All BSpline<6> grids now share ONE coarse
// "BSpline<6>" Cache4 (RadialType()==Name(), no grid in the key), so correctness rests on the lossless
// SplineGrouper (full knot-vector key) plus per-grid GridData routing.  Two DIFFERENT grids (same N so the
// ERI4 blocks are comparable, different rmin/rmax) must therefore yield DISTINCT integrals -- if the
// grouper aliased their splines, or Create() routed to the wrong grid's GL/Rk tables, JB would collapse
// onto JA's cached values and fnorm(JA,JB) would go to 0.
TEST_F(Cache4Tests,HF2_BS_GridKeyed)
{
    delete bs1;
    nlohmann::json jsA={{"N", 5}, {"rmin", 0.25}, {"rmax", 4.0}}; jsA["type"]=Type::BSpline6;
    nlohmann::json jsB={{"N", 5}, {"rmin", 0.20}, {"rmax", 3.0}}; jsB["type"]=Type::BSpline6;
    BasisSet::Real_BS* bsA=Factory(jsA,2);
    BasisSet::Real_BS* bsB=Factory(jsB,2);   // registers into the SAME "BSpline<6>" Cache4 as bsA
    using BasisSet::Real_ERI4_OIBS;
    auto a=bsA->Iterate<Real_ERI4_OIBS>().begin();
    auto b=bsB->Iterate<Real_ERI4_OIBS>().begin();
    const ERI4& JA=(*a)->Direct(**a);
    const ERI4& JB=(*b)->Direct(**b);
    EXPECT_GT(fnorm(JA,JB),1e-3);            // distinct grids -> distinct integrals (no cross-grid aliasing)
    delete bsA;
    delete bsB;
}

// Each atom basis owns its Rk table (S5a).  (The exponent Pool, src/BasisSet/Radial/Factory.C, still makes a lighter element's exponents
// a bit-identical subset of a heavier one's; the process-wide Cache4 that turned that into cross-element reuse is gone with S5a.)
TEST_F(Cache4Tests,HF2_SL_OneTablePerAtom)
{
    TestOneTablePerAtom(Type::Slater,"SL",92,68);
}
TEST_F(Cache4Tests,HF2_SG_OneTablePerAtom)
{
    TestOneTablePerAtom(Type::Gaussian,"SG",54,36);
}


// TEST_F(Cache4Tests,Caching)
// {
//     auto cache=BasisSet::theGlobalCache;
//     // EXPECT_NE(cache,NULL);
//     // assert(cache);

//     auto s=new Evaluator({1,2,4,8},0);
//     auto p=new Evaluator({.5,1,2.0,4},1);
//     auto d=new Evaluator({.25,.5,1.0,2},2);
//     auto f=new Evaluator({.25,.5,1.0},3);

//     cache->Register(s);
//     cache->Register(p);
//     cache->Register(d);
//     cache->Register(f);

//     const Cache4* sc=cache->GetCache4(s->RadialType());
//     EXPECT_EQ(sc,cache->GetCache4(p->RadialType()));
//     EXPECT_EQ(sc,cache->GetCache4(d->RadialType()));
//     EXPECT_EQ(sc,cache->GetCache4(f->RadialType()));

//     const Gaussian_Cache4* sgc=dynamic_cast<const Gaussian_Cache4*>(sc);
//     EXPECT_TRUE(sgc);

//     const ExponentGrouper& gr=GetGrouper(sgc);

//     EXPECT_THAT(gr.unique_esv,::testing::ElementsAre(1,2,4,8,.5,.25));
//     // std::cout << "Exponents: " << gr.unique_esv << std::endl;
//     // std::cout << "maxls: " << maxls(sgc) << std::endl;
//     // std::cout << "s->es_indices: " << es_indices(s) << std::endl;
//     EXPECT_THAT(gr.unique_esv,::testing::ElementsAre(1,2,4,8,.5,.25));
//     EXPECT_THAT(maxls(sgc),::testing::ElementsAre(3,2,1,0,3,3));
//     EXPECT_THAT(es_indices(s),::testing::ElementsAre(0,1,2,3));
//     EXPECT_THAT(es_indices(p),::testing::ElementsAre(4,0,1,2));
//     EXPECT_THAT(es_indices(d),::testing::ElementsAre(5,4,0,1));
//     EXPECT_THAT(es_indices(f),::testing::ElementsAre(5,4,0));

// }

size_t ClosedShellZs[]={2,10,18,36,46,63,70,80,88};
size_t   OpenShellZs[]={5,6,7,8,9,11,13,14,15,16,17,19,21,22,23,24,25,26,27,28,29,39,40,41,42,43,58,64,91,92};

TEST_F(Cache4Tests,HF2_SG_Direct_ClosedShell)
{
    for (auto Z:ClosedShellZs)
    {
        InitGaussian(Z);
        TestDirect(2.2e-15);
    }
}
TEST_F(Cache4Tests,HF2_SG_Exchange_ClosedShell)
{
   for (auto Z:ClosedShellZs)
    {
        InitGaussian(Z);
        TestExchange(2.2e-15);
    }
}
TEST_F(Cache4Tests,HF2_SG_Direct_OpenShell)
{
    for (auto Z:OpenShellZs)
    {
        InitGaussian(Z);
        TestDirect(2.2e-15);
    }
}
TEST_F(Cache4Tests,HF2_SG_Exchange_OpenShell)
{
    for (auto Z:OpenShellZs)
    {
        InitGaussian(Z);
        TestExchange(2.2e-15);
    }
}
TEST_F(Cache4Tests,HF2_SL_Direct_ClosedShell)
{
    for (auto Z:ClosedShellZs)
    {
        InitSlater(Z);
        TestDirect(2.2e-15);
    }
}
TEST_F(Cache4Tests,HF2_SL_Exchange_ClosedShell)
{
    for (auto Z:ClosedShellZs)
    {
        InitSlater(Z);
        TestExchange(2.2e-15);
    }
}
TEST_F(Cache4Tests,HF2_SL_Direct_OpenShell)
{
    for (auto Z:OpenShellZs)
    {
        InitSlater(Z);
        TestDirect(2.2e-15);
    }
}
TEST_F(Cache4Tests,HF2_SL_Exchange_OpenShell)
{
     for (auto Z:OpenShellZs)
    {
        InitSlater(Z);
        TestExchange(2.2e-15);
    }
}

TEST_F(Cache4Tests,HF2_BS_Direct_ClosedShell)
{
    for (auto Z:ClosedShellZs)
    {
        InitBSpline6(Z);
        TestDirect(2.2e-15,false);
    }
}
TEST_F(Cache4Tests,HF2_BS_Exchange_ClosedShell)
{
    for (auto Z:ClosedShellZs)
    {
        InitBSpline6(Z);
        TestExchange(2.2e-15,false);
    }
}
TEST_F(Cache4Tests,HF2_BS_Direct_OpenShell)
{
    for (auto Z:OpenShellZs)
    {
        InitBSpline6(Z);
        TestDirect(2.2e-15,false);
    }
}
TEST_F(Cache4Tests,HF2_BS_Exchange_OpenShell)
{
    for (auto Z:OpenShellZs)
    {
        InitBSpline6(Z);
        TestExchange(2.2e-15,false);
    }
}

