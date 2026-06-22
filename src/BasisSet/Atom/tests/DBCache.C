// File: UnitTests/BasisSet_Atom1.C  Unit test the Atom IBS Evaluators
#include "gtest/gtest.h"
#include <iostream>
#include <nlohmann/json.hpp>
#include <blaze/math/expressions/DMatDMatEqualExpr.h> //op== inside gtest header.

using std::cout;
using std::endl;

import qchem.BasisSet.Atom.Factory;
import qchem.BasisSet.Orbital_HF_IBS;
import qchem.BasisSet.Orbital_DFT_IBS;
import qchem.Blaze;

using namespace BasisSet::Atom;
using enum Type;
using enum BasisSetAccuracy;

class DBCach1Tests : public ::testing::Test
{
public:
    DBCach1Tests() 
        : N(3),Z(1)
        , cl_hydrogen    (new Atom(1,0.0,Vector3D(0,0,0)))
        , cl_hydrogen_100(new Atom(1,0.0,Vector3D(1,0,0)))
        , cl_helium      (new Atom(2,0.0,Vector3D(0,0,0)))
    {
    }
    ~DBCach1Tests()
    {
        delete cl_hydrogen;
        delete cl_hydrogen_100;
        delete cl_helium;
        delete bs1;
        delete bs2;
    }
    void Init(BasisSet::Atom::Type type)
    {
        bs1=Factory(Low,type,Z);
        bs2=Factory(Low,type,Z);
    }
  
    void TestOverlap() const
    {
        using OIBS=BasisSet::Real_OIBS;
        auto ibs2=bs2->Iterate<OIBS>().begin();
        for (auto ibs1:bs1->Iterate<OIBS>())
        {
            auto& S1=ibs1->Overlap();
            auto& S2=(*ibs2)->Overlap();
            EXPECT_EQ(S1,S2);
            EXPECT_EQ(&S1,&S2);
            ++ibs2;
        }
    }
    void TestKinetic() const
    {
        using OIBS=BasisSet::Real_OIBS;
        auto ibs2=bs2->Iterate<OIBS>().begin();
        for (auto ibs1:bs1->Iterate<OIBS>())
        {
            // Kinetic() here is just a representative cached 1e matrix (the <p^2> block) used to
            // check the cache returns the SAME object for the same key -- its physics is incidental.
            auto& S1=ibs1->Kinetic();
            auto& S2=(*ibs2)->Kinetic();
            EXPECT_EQ(S1,S2);
            EXPECT_EQ(&S1,&S2);
            ++ibs2;
        }
    }
    void TestNuclear() const
    {
        using OIBS=BasisSet::Real_OIBS;
        auto ibs2=bs2->Iterate<OIBS>().begin();
        for (auto ibs1:bs1->Iterate<OIBS>())
        {
            auto& S1=ibs1->Nuclear(cl_hydrogen);
            auto& S2=(*ibs2)->Nuclear(cl_hydrogen);
            auto& S3=(*ibs2)->Nuclear(cl_hydrogen_100);
            auto& S4=(*ibs2)->Nuclear(cl_helium);
            EXPECT_EQ(S1,S2);
            EXPECT_EQ(S1,S3);
            EXPECT_NE(S1,S4);
            EXPECT_EQ(&S1,&S2);
            EXPECT_NE(&S1,&S3);
            EXPECT_NE(&S1,&S4);
            ++ibs2;
        }
    }
    void TestOverlap3C(double eps) const
    {
        using BasisSet::Real_DFT_OIBS;
        auto ibs2=bs2->Iterate<Real_DFT_OIBS>().begin();
        for (auto ibs1:bs1->Iterate<Real_DFT_OIBS>())
        {
            auto ff=ibs1->CreateCDFitBasisSet(cl_hydrogen);
            const ERI3<double>& E1=ibs1->Overlap3C(*ff);
            const ERI3<double>& E2=(*ibs2)->Overlap3C(*ff);
            EXPECT_EQ(E1,E2);
            EXPECT_EQ(&E1,&E2);
            ++ibs2;
        }

    }
    void TestRepulsion3C(double eps) const
        {
            using BasisSet::Real_DFT_OIBS;
            auto ibs2=bs2->Iterate<Real_DFT_OIBS>().begin();
            for (auto ibs1:bs1->Iterate<Real_DFT_OIBS>())
            {
                auto ff=ibs1->CreateCDFitBasisSet(cl_hydrogen);
                const ERI3<double>& E1=ibs1->Repulsion3C(*ff);
                const ERI3<double>& E2=(*ibs2)->Repulsion3C(*ff);
                EXPECT_EQ(E1,E2);
                EXPECT_EQ(&E1,&E2);
                ++ibs2;
            }

        }
    void TestDirect(double eps) const
        {
            cout << *bs1 << endl;
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
                    EXPECT_NEAR(fnorm(J1,J1ba.Transpose()),0.0,eps);
                    EXPECT_NEAR(fnorm(J2,J2ba.Transpose()),0.0,eps);
                    ++ibs22;
                    // cout << std::setprecision(12) << "J1(0,0)(0,1)=" << J1(0,0)(0,1) << std::endl;
                    // cout << std::setprecision(12) << "J1(0,1)(0,0)=" << J1(0,1)(0,0) << std::endl;
                }
                ++ibs21;
            }

        }
    void TestExchange(double eps) const
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
                EXPECT_NEAR(fnorm(K1,K1ba.Transpose()),0.0,eps);
                EXPECT_NEAR(fnorm(K2,K2ba.Transpose()),0.0,eps);
                ++ibs22;
            }
            ++ibs21;
        }
    }
    size_t N,Z;
    Structure *cl_hydrogen,*cl_hydrogen_100,*cl_helium;
    BasisSet::Real_BS *bs1,*bs2;
};




TEST_F(DBCach1Tests,BSplineOverlap)
{
    Init(BSpline6);
    TestOverlap();
}
TEST_F(DBCach1Tests,GaussianOverlap)
{
    Init(Gaussian);
    TestOverlap();
}
TEST_F(DBCach1Tests,SlaterOverlap)
{
    Init(Slater);
    TestOverlap();
}

TEST_F(DBCach1Tests,BSplineKinetic)
{
    Init(BSpline6);
    TestKinetic();
}
TEST_F(DBCach1Tests,GaussianKinetic)
{
    Init(Gaussian);
    TestKinetic();
}
TEST_F(DBCach1Tests,SlaterKinetic)
{
    Init(Slater);
    TestKinetic();
}
TEST_F(DBCach1Tests,BSplineNuclear)
{
    Init(BSpline6);
    TestNuclear();
}
TEST_F(DBCach1Tests,GaussianNuclear)
{
    Init(Gaussian);
    TestNuclear();
}
TEST_F(DBCach1Tests,SlaterNuclear)
{
    Init(Slater);
    TestNuclear();
}
TEST_F(DBCach1Tests,GaussianOverlap3C)
{
    Init(Gaussian);
    TestOverlap3C(1e-14);
}
TEST_F(DBCach1Tests,SlaterOverlap3C)
{
    Init(Slater);
    TestOverlap3C(1e-14);
}
TEST_F(DBCach1Tests,GaussianRepulsion3C)
{
    Init(Gaussian);
    TestRepulsion3C(1e-14);
}
TEST_F(DBCach1Tests,SlaterRepulsion3C)
{
    Init(Slater);
    TestRepulsion3C(1e-14);
}


//
//  For efficiency legacy calculations do Jcdab=Jabcd.Transpose().   
//  The BasisSet version recalculates Jcdab directly and does not yet do
//  the Transpose call.  For Guassian and Slater this results in a minor error <1e-15.  For BSpline near the origin
//  we get much bigger errors.  Perhaps this related to problems encountered when trying to reproduce Saitos results with
//  rmin=0.0001.
//
TEST_F(DBCach1Tests,BSplineDirect)
{
    Init(BSpline6);
    TestDirect(6e-13);
}
TEST_F(DBCach1Tests,BSplinerDirect)
{
    Init(BSpliner6);
    TestDirect(1.2e-13);
}
TEST_F(DBCach1Tests,GaussianDirect)
{
    Init(Gaussian);
    TestDirect(7e-14);
}
TEST_F(DBCach1Tests,SlaterDirect)
{
    Init(Slater);
    TestDirect(4e-15);
}
TEST_F(DBCach1Tests,BSplineExchange)
{
    Init(BSpline6);
    TestExchange(1.4e-13);
}
TEST_F(DBCach1Tests,BSplinerExchange)
{
    Init(BSpliner6);
    TestExchange(6e-14);
}
TEST_F(DBCach1Tests,GaussianExchange)
{
    Init(Gaussian);
    TestExchange(6e-14);
}
TEST_F(DBCach1Tests,SlaterExchange)
{
    Init(Slater);
    TestExchange(4e-15);
}
