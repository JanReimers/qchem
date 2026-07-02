// File: UnitTests/A_Cache4.C  Unit test the Atom Cach4 system for storing 4 index radial Slater integrals.
#include "gtest/gtest.h"
#include <gmock/gmock.h>
#include <iostream>
#include <map>
#include <string>
#include <memory>
#include "../../src/forward.H"
#include <nlohmann/json.hpp>
using std::cout;
using std::endl;

import qchem.BasisSet.Atom.Evaluators; 
import qchem.BasisSet.Atom.Evaluators.Gaussian.IBS; 
import qchem.BasisSet.Atom.Evaluators.Internal.Rk;
import qchem.BasisSet.Atom.Evaluators.Internal.Grouper;
import qchem.BasisSet.Atom.Evaluators.Gaussian.Internal.Rk; 
import qchem.BasisSet.Atom.Evaluators.Internal.ExponentialEvaluator;
import qchem.stl_io;
import qchem.BasisSet.Atom.Factory;
import qchem.BasisSet.Orbital_HF_IBS;
import qchem.BasisSet.Internal.DB_Cache_RAM;   // theCache<double>() + GetCache4(RadialType)
import qchem.BasisSet.Internal.Cache4;         // Cache4::Lookups()/Inserts()
using namespace qchem;

using namespace qchem::BasisSet::Atom;
using namespace qchem::BasisSet::Atom::Evaluators;
using enum BasisSetAccuracy;
using GCache4=BasisSet::Atom::Evaluators::Gaussian::Gaussian_Cache4;

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
        bs1=BasisSet::Atom::Factory(N3,type1,Z);
        bs2=BasisSet::Atom::Factory(N3,type2,Z);
    }
    void InitBSpline6(size_t Z) {Init(Z,Type::BSpline6,Type::BSpline6);}
    void InitGaussian(size_t Z) {Init(Z,Type::Gaussian,Type::Gaussian);}
    void InitSlater  (size_t Z) {Init(Z,Type::Slater  ,Type::Slater);}

    void TestDirect(double eps, bool testTranspose=true)
    {
        using BasisSet::Real_HF_OIBS;
        auto ibs21=bs2->Iterate<Real_HF_OIBS>().begin();
        for (auto ibs11:bs1->Iterate<Real_HF_OIBS>())
        {
            auto ibs22=bs2->Iterate<Real_HF_OIBS>().begin();
            for (auto ibs12:bs1->Iterate<Real_HF_OIBS>())
            {
                const ERI4& J1=ibs11->Direct(*ibs12);
                const ERI4& J2=(*ibs21)->Direct(**ibs22);
                EXPECT_EQ(J1,J2);
                EXPECT_EQ(&J1,&J2);
                const ERI4& J1ba=ibs12->Direct(*ibs11);
                const ERI4& J2ba=(*ibs22)->Direct(**ibs21);
                if (testTranspose)
                {
                    EXPECT_NEAR(fnorm(J1,J1ba.Transpose()),0.0,eps);
                    EXPECT_NEAR(fnorm(J2,J2ba.Transpose()),0.0,eps);
                }
                ++ibs22;
            }
            ++ibs21;
        }

    }
    void TestExchange(double eps, bool testTranspose=true)
    {
        using BasisSet::Real_HF_OIBS;
        auto ibs21=bs2->Iterate<Real_HF_OIBS>().begin();
        for (auto ibs11:bs1->Iterate<Real_HF_OIBS>())
        {
            auto ibs22=bs2->Iterate<Real_HF_OIBS>().begin();
            for (auto ibs12:bs1->Iterate<Real_HF_OIBS>())
            {
                const ERI4& K1=ibs11->Exchange(*ibs12);
                const ERI4& K2=(*ibs21)->Exchange(**ibs22);
                EXPECT_EQ(K1,K2);
                EXPECT_EQ(&K1,&K2);
                const ERI4& K1ba=ibs12->Exchange(*ibs11);
                const ERI4& K2ba=(*ibs22)->Exchange(**ibs21);
                if (testTranspose)
                {
                    EXPECT_NEAR(fnorm(K1,K1ba.Transpose()),0.0,eps);
                    EXPECT_NEAR(fnorm(K2,K2ba.Transpose()),0.0,eps);
                }
                ++ibs22;
            }
            ++ibs21;
        }

    }
    
    void TestDirect(size_t Z,nlohmann::json js)
    {
        bs1=Factory(js,Z);
        using BasisSet::Real_HF_OIBS;
        for (auto ibs:bs1->Iterate<Real_HF_OIBS>())
        {
            // std::cout << *ibs;
            ibs->Direct(*ibs);
        }
        delete bs1;
    }
 
    

    // Trigger every 4-index Rk Create() for a basis by evaluating Direct+Exchange over all shell pairs,
    // then dispose of the basis.  The Rk integrals it computed persist in the process-wide Cache4.
    void ExerciseAllRk(BasisSet::Real_BS* bs)
    {
        using BasisSet::Real_HF_OIBS;
        for (auto i:bs->Iterate<Real_HF_OIBS>())
            for (auto j:bs->Iterate<Real_HF_OIBS>())
            {
                i->Direct  (*j);
                i->Exchange(*j);
            }
        delete bs;
    }

    // The exponent-Pool payoff: a heavy atom warms the shared Rk cache for every lighter atom of the
    // same (family, accuracy) because the lighter atom's exponents are a bit-identical SUBSET.  Warm
    // with the heavy atom, snapshot the cache, then run the light atom and confirm ~100% reuse (a few
    // low-l-only Rks are re-created because the light atom's higher-l shells over-evict them on Register
    // -- harmless and recreated correctly -- so we bound the miss rate rather than demand exactly zero).
    void TestCrossElementReuse(Type type, const std::string& radialType, size_t Zheavy, size_t Zlight)
    {
        auto& cache=BasisSet::theCache<double>();
        ExerciseAllRk(Factory(Low,type,Zheavy));            // heavy atom populates the shared Rk cache
        const Cache4* c4=cache.GetCache4(radialType);
        ASSERT_TRUE(c4);
        const size_t inserts0=c4->Inserts(), lookups0=c4->Lookups();
        ExerciseAllRk(Factory(Low,type,Zlight));            // light atom: a strict subset of the heavy exponents
        const size_t dInserts=c4->Inserts()-inserts0;
        const size_t dLookups=c4->Lookups()-lookups0;
        const double reuse=dLookups ? 100.0*(1.0-double(dInserts)/double(dLookups)) : 0.0;
        std::cout << "[pool] " << radialType << " Z" << Zlight << "-after-Z" << Zheavy
                  << ": dLookups=" << dLookups << " dInserts=" << dInserts
                  << " reuse=" << reuse << "%" << std::endl;
        EXPECT_GT(dLookups,10000u);                          // light atom really did substantial Rk work
        EXPECT_GT(reuse,99.9);                               // ...and nearly all of it hit the heavy atom's cache
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
    using BasisSet::Real_HF_OIBS;
    auto a=bsA->Iterate<Real_HF_OIBS>().begin();
    auto b=bsB->Iterate<Real_HF_OIBS>().begin();
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
    using BasisSet::Real_HF_OIBS;
    auto a=bsA->Iterate<Real_HF_OIBS>().begin();
    auto b=bsB->Iterate<Real_HF_OIBS>().begin();
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
    using BasisSet::Real_HF_OIBS;
    auto a=bsA->Iterate<Real_HF_OIBS>().begin();
    auto b=bsB->Iterate<Real_HF_OIBS>().begin();
    const ERI4& JA=(*a)->Direct(**a);
    const ERI4& JB=(*b)->Direct(**b);
    EXPECT_GT(fnorm(JA,JB),1e-3);            // distinct grids -> distinct integrals (no cross-grid aliasing)
    delete bsA;
    delete bsB;
}

// Cross-element Rk reuse (the "exponent Pool" payoff, src/BasisSet/Atom/Factory.C).
// SlaterExponents()/GaussianExponents() build ONE universal even-tempered pool: emin/beta/emax are fixed
// constants, so e*=beta yields a bit-identical geometric sequence for EVERY Z; only the COUNT is trimmed
// by Z (largest exponents dropped; High keeps the smallest NZ, Low strides that prefix).  Hence a lighter
// element's exponent set is a strict, bit-identical subset of a heavier one's.  So once Uranium(92) has
// filled the shared "SL"/"SG" Rk Cache4, Erbium(68) reuses every entry -- zero new inserts, 100% reuse.
// (Medium is deliberately NOT used here: it trims small exponents proportionally, so the lighter atom
// keeps a small exponent the heavier one dropped and the subset relation breaks.)
TEST_F(Cache4Tests,HF2_SL_CrossElementReuse)
{
    TestCrossElementReuse(Type::Slater,"SL",92,68);         // Uranium warms the pool for Erbium
}
TEST_F(Cache4Tests,HF2_SG_CrossElementReuse)
{
    TestCrossElementReuse(Type::Gaussian,"SG",54,36);       // lighter pair keeps the (denser) Gaussian run cheap
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

