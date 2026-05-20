// File: UnitTests/BasisSet_DHF.C  Unit test the Atom RKB basis sets
#include "gtest/gtest.h"
#include <iostream>
// #include <cmath>
#include <blaze/Math.h>
#include <nlohmann/json.hpp>
using std::cout;
using std::endl;

import qchem.BasisSet.DB_Cache;
import qchem.BasisSet.Atom.Factory;
import qchem.BasisSet.Orbital_DHF_IBS;


class Basis1_RKB_Tests : public ::testing::Test
{
public:
    Basis1_RKB_Tests() 
        : N(3),Z(1)
        , cl_hydrogen    (new Atom(1,0.0,Vector3D(0,0,0)))
        , cl_hydrogen_100(new Atom(1,0.0,Vector3D(1,0,0)))
        , cl_helium      (new Atom(2,0.0,Vector3D(0,0,0)))
    {
        // BasisSet::theGlobalCache=new BasisSet::IntegralsCache_RAM<double>();       
    }
    ~Basis1_RKB_Tests()
    {
        delete cl_hydrogen;
        delete cl_hydrogen_100;
        delete cl_helium;
        delete bs1;
        delete bs2;
        // delete BasisSet::theGlobalCache;
    }
    void Init(nlohmann::json js)
    {
        bs1=BasisSet::Atom::Factory(js,Z);
        bs2=BasisSet::Atom::Factory(js,Z);
    }
    // void InitBSpline6() {Init({{"type",BasisSet::Atom::Type::BSpline6},{"N", N}, {"rmin", 0.1}, {"rmax", 10}},BasisSetAtom::Type::BSpline6);}
    void InitGaussian() {Init({{"type",BasisSet::Atom::Type::Gaussian_RKB},{"N", N}, {"emin", 0.1}, {"emax", 10}});}
    void InitSlater  () {Init({{"type",BasisSet::Atom::Type::Slater_RKB  },{"N", N}, {"emin", 0.1}, {"emax", 10}});}
    using OIBS=BasisSet::Real_ORKB;
    void TestOverlap() const
    {
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
    void TestRestMass() const
    {
        auto ibs2=bs2->Iterate<OIBS>().begin();
        for (auto ibs1:bs1->Iterate<OIBS>())
        {
            auto& S1=ibs1->RestMass();
            auto& S2=(*ibs2)->RestMass();
            EXPECT_EQ(S1,S2);
            EXPECT_EQ(&S1,&S2);
            ++ibs2;
        }
    }
    
    size_t N,Z;
    Cluster *cl_hydrogen,*cl_hydrogen_100,*cl_helium;
    BasisSet::Real_BS *bs1,*bs2;
};




TEST_F(Basis1_RKB_Tests,GaussianOverlap)
{
    InitGaussian();
    TestOverlap();
}
TEST_F(Basis1_RKB_Tests,GaussianKinetic)
{
    InitGaussian();
    TestKinetic();
}
TEST_F(Basis1_RKB_Tests,GaussianNuclear)
{
    InitGaussian();
    TestNuclear();
}
TEST_F(Basis1_RKB_Tests,GaussianRestMass)
{
    InitGaussian();
    TestRestMass();
}
TEST_F(Basis1_RKB_Tests,SlaterOverlap)
{
    InitSlater();
    TestOverlap();
}
TEST_F(Basis1_RKB_Tests,SlaterKinetic)
{
    InitSlater();
    TestKinetic();
}
TEST_F(Basis1_RKB_Tests,SlaterNuclear)
{
    InitSlater();
    TestNuclear();
}
TEST_F(Basis1_RKB_Tests,SlaterRestMass)
{
    InitSlater();
    TestRestMass();
}
