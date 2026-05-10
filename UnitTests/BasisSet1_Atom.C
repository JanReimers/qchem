// File: UnitTests/BasisSet_Atom1.C  Unit test the Atom IBS Evaluators
#include "gtest/gtest.h"
#include <iostream>
// #include <cmath>
#include <blaze/Math.h>
#include <nlohmann/json.hpp>
using std::cout;
using std::endl;

import qchem.BasisSet1.DB_Cache;
import qchem.BasisSet1.Atom.Factory;
import qchem.BasisSet1.Orbital_HF_IBS;
import qchem.BasisSet1.Orbital_DFT_IBS;

// Legacy BS imports
// import qchem.BasisSet1.Atom.Factory;
// import qchem.Orbital_HF_IBS;
// import qchem.Orbital_DFT_IBS;


class DBCach1Tests : public ::testing::Test
{
public:
    DBCach1Tests() 
        : N(3),Z(86)
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
    void Init(nlohmann::json js)
    {
        bs1=BasisSet1::Atom::Factory(js,Z);
        bs2=BasisSet1::Atom::Factory(js,Z);
    }
    void InitBSpline6() {Init({{"type",BasisSet1::Atom::Type::BSpline6},{"N", N}, {"rmin", 0.1}, {"rmax", 10}});}
    void InitGaussian() {Init({{"type",BasisSet1::Atom::Type::Gaussian},{"N", N}, {"emin", 0.1}, {"emax", 10}});}
    void InitSlater  () {Init({{"type",BasisSet1::Atom::Type::Slater  },{"N", N}, {"emin", 0.1}, {"emax", 10}});}

    void TestOverlap() const
    {
        using OIBS=BasisSet1::Real_OIBS;
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
        using OIBS=BasisSet1::Real_OIBS;
        auto ibs2=bs2->Iterate<OIBS>().begin();
        for (auto ibs1:bs1->Iterate<OIBS>())
        {
            auto& S1=ibs1->Kinetic();
            auto& S2=(*ibs2)->Kinetic();
            EXPECT_EQ(S1,S2);
            EXPECT_EQ(&S1,&S2);
            ++ibs2;
        }
    }
    void TestNuclear() const
    {
        using OIBS=BasisSet1::Real_OIBS;
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
        using BasisSet1::Real_DFT_OIBS;
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
            using BasisSet1::Real_DFT_OIBS;
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
                    ++ibs22;
                }
                ++ibs21;
            }

        }
    void TestExchange(double eps) const
    {
        
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
                ++ibs22;
            }
            ++ibs21;
        }
    }
    size_t N,Z;
    Cluster *cl_hydrogen,*cl_hydrogen_100,*cl_helium;
    BasisSet1::Real_BS *bs1,*bs2;
};




TEST_F(DBCach1Tests,BSplineOverlap)
{
    InitBSpline6();
    TestOverlap();
}
TEST_F(DBCach1Tests,GaussianOverlap)
{
    InitGaussian();
    TestOverlap();
}
TEST_F(DBCach1Tests,SlaterOverlap)
{
    InitSlater();
    TestOverlap();
}

TEST_F(DBCach1Tests,BSplineKinetic)
{
    InitBSpline6();
    TestKinetic();
}
TEST_F(DBCach1Tests,GaussianKinetic)
{
    InitGaussian();
    TestKinetic();
}
TEST_F(DBCach1Tests,SlaterKinetic)
{
    InitSlater();
    TestKinetic();
}
TEST_F(DBCach1Tests,BSplineNuclear)
{
    InitBSpline6();
    TestNuclear();
}
TEST_F(DBCach1Tests,GaussianNuclear)
{
    InitGaussian();
    TestNuclear();
}
TEST_F(DBCach1Tests,SlaterNuclear)
{
    InitSlater();
    TestNuclear();
}
TEST_F(DBCach1Tests,GaussianOverlap3C)
{
    InitGaussian();
    TestOverlap3C(1e-14);
}
TEST_F(DBCach1Tests,SlaterOverlap3C)
{
    InitSlater();
    TestOverlap3C(1e-14);
}
TEST_F(DBCach1Tests,GaussianRepulsion3C)
{
    InitGaussian();
    TestRepulsion3C(1e-14);
}
TEST_F(DBCach1Tests,SlaterRepulsion3C)
{
    InitSlater();
    TestRepulsion3C(1e-14);
}


//
//  For efficiency legacy calculations do Jcdab=Jabcd.Transpose().   
//  The BasisSet1 version recalculates Jcdab directly and does not yet do
//  the Transpose call.  For Guassian and Slater this results in a minor error <1e-15.  For BSpline near the origin
//  we get much bigger errors.  Perhaps this related to problems encountered when trying to reproduce Saitos results with
//  rmin=0.0001.
//
TEST_F(DBCach1Tests,BSplineDirect)
{
    InitBSpline6();
    TestDirect(1e-14);
}
TEST_F(DBCach1Tests,GaussianDirect)
{
    InitGaussian();
    TestDirect(3e-15);
}
TEST_F(DBCach1Tests,SlaterDirect)
{
    InitSlater();
    TestDirect(4e-15);
}
TEST_F(DBCach1Tests,BSplineExchange)
{
    InitBSpline6();
    TestExchange(1e-14);
}
TEST_F(DBCach1Tests,GaussianExchange)
{
    InitGaussian();
    TestExchange(4e-15);
}
TEST_F(DBCach1Tests,SlaterExchange)
{
    InitSlater();
    TestExchange(4e-15);
}
