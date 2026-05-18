// File: UnitTests/A_Cache4.C  Unit test the Atom Cach4 system for storing 4 index redial Slater integrals.
#include "gtest/gtest.h"
#include <gmock/gmock.h>
#include <iostream>
#include <map>
#include <string>
#include <memory>
#include "../src/forward.H"
// #include <blaze/Math.h>
#include <nlohmann/json.hpp>
using std::cout;
using std::endl;

import qchem.BasisSet1.DB_Cache;
import qchem.BasisSet1.Atom.Evaluators.IBS; 
import qchem.BasisSet1.Atom.Evaluators.Gaussian.IBS; 
import qchem.BasisSet1.Atom.Evaluators.Internal.Rk;
import qchem.BasisSet1.Atom.Evaluators.Internal.ExponentGrouper;
import qchem.BasisSet1.Atom.Evaluators.Gaussian.Internal.Rk; 
import qchem.stl_io;
// import qchem.BasisSet1.Internal.Cache4;
import qchem.BasisSet1.Atom.Factory;
// import qchem.BasisSet1.Atom.BasisSet;
import qchem.Symmetry.AtomEC;
import qchem.BasisSet1.Orbital_HF_IBS;

class Cache4Tests : public ::testing::Test
{
public:
    Cache4Tests() 
    {
        if (BasisSet1::theGlobalCache==0)
            BasisSet1::theGlobalCache=new BasisSet1::IntegralsCache_RAM<double>(true);     

    }

    const ExponentGrouper& GetGrouper(const Gaussian_Cache4* gc)
    {
        assert(gc);
        return gc->grouper; //friend access.
    }

    const std::vector<size_t>& maxls(const Gaussian_Cache4* gc)
    {
        return GetGrouper(gc).maxls;
    }

    const std::vector<size_t>& es_indices(const IBS_Evaluator* e)
    {
        return e->es_indices;
    }

};


TEST_F(Cache4Tests,Caching)
{
    auto cache=BasisSet1::theGlobalCache;
    // EXPECT_NE(cache,NULL);
    // assert(cache);

    auto s=new Gaussian_IBS_Evaluator({1,2,4,8},0);
    auto p=new Gaussian_IBS_Evaluator({.5,1,2.0,4},1);
    auto d=new Gaussian_IBS_Evaluator({.25,.5,1.0,2},2);
    auto f=new Gaussian_IBS_Evaluator({.25,.5,1.0},3);

    cache->Register(s);
    cache->Register(p);
    cache->Register(d);
    cache->Register(f);

    const Cache41* sc=cache->GetCache4(s->RadialType());
    EXPECT_EQ(sc,cache->GetCache4(p->RadialType()));
    EXPECT_EQ(sc,cache->GetCache4(d->RadialType()));
    EXPECT_EQ(sc,cache->GetCache4(f->RadialType()));

    const Gaussian_Cache4* sgc=dynamic_cast<const Gaussian_Cache4*>(sc);
    EXPECT_TRUE(sgc);

    const ExponentGrouper& gr=GetGrouper(sgc);

    EXPECT_THAT(gr.unique_esv,::testing::ElementsAre(1,2,4,8,.5,.25));
    // std::cout << "Exponents: " << gr.unique_esv << std::endl;
    // std::cout << "maxls: " << maxls(sgc) << std::endl;
    // std::cout << "s->es_indices: " << es_indices(s) << std::endl;
    EXPECT_THAT(gr.unique_esv,::testing::ElementsAre(1,2,4,8,.5,.25));
    EXPECT_THAT(maxls(sgc),::testing::ElementsAre(3,2,1,0,3,3));
    EXPECT_THAT(es_indices(s),::testing::ElementsAre(0,1,2,3));
    EXPECT_THAT(es_indices(p),::testing::ElementsAre(4,0,1,2));
    EXPECT_THAT(es_indices(d),::testing::ElementsAre(5,4,0,1));
    EXPECT_THAT(es_indices(f),::testing::ElementsAre(5,4,0));

}

TEST_F(Cache4Tests,BasisSet_HF2_Direct)
{
    double eps=2e-15;
    Atom_EC ec(2);
    nlohmann::json js1={{"type",BasisSet1::Atom::Type::Gaussian },{"N", 3}, {"emin", 0.5}, {"emax", 2}};
    nlohmann::json js2={{"type",BasisSet1::Atom::Type::Gaussian2},{"N", 3}, {"emin", 0.5}, {"emax", 2}};
    auto bs1=BasisSet1::Atom::Factory(js1,ec);
    auto bs2=BasisSet1::Atom::Factory(js2,ec);
    using BasisSet1::Real_HF_OIBS;
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
            EXPECT_NEAR(fnorm(J1,J1ba.Transpose()),0.0,eps);
            EXPECT_NEAR(fnorm(J2,J2ba.Transpose()),0.0,eps);
            ++ibs22;
        }
        ++ibs21;
    }
}

TEST_F(Cache4Tests,BasisSet_HF2_Exchange)
{
    double eps=2e-15;
    Atom_EC ec(2);
    nlohmann::json js1={{"type",BasisSet1::Atom::Type::Gaussian },{"N", 3}, {"emin", 0.5}, {"emax", 2}};
    nlohmann::json js2={{"type",BasisSet1::Atom::Type::Gaussian2},{"N", 3}, {"emin", 0.5}, {"emax", 2}};
    auto bs1=BasisSet1::Atom::Factory(js1,ec);
    auto bs2=BasisSet1::Atom::Factory(js2,ec);
    using BasisSet1::Real_HF_OIBS;
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
            EXPECT_NEAR(fnorm(K1,K1ba.Transpose()),0.0,eps);
            EXPECT_NEAR(fnorm(K2,K2ba.Transpose()),0.0,eps);
            ++ibs22;
        }
        ++ibs21;
    }
}