// File: UnitTests/BasisSet_DHF.C  Unit test the Atom RKB basis sets
#include "gtest/gtest.h"
#include <iostream>
// #include <cmath>
#include <blaze/Math.h>
#include <nlohmann/json.hpp>
using std::cout;
using std::endl;

import qchem.BasisSet1.DB_Cache;
import qchem.BasisSet1.Atom.Factory;
import qchem.BasisSet1.Orbital_DHF_IBS;

// Legacy BS imports
import qchem.BasisSet.Atom.Factory;
import qchem.Orbital_DHF_IBS;


class Basis1_RKB_Tests : public ::testing::Test
{
public:
    Basis1_RKB_Tests() 
        : N(3),Z(1)
        , cl_hydrogen    (new Atom(1,0.0,Vector3D(0,0,0)))
        , cl_hydrogen_100(new Atom(1,0.0,Vector3D(1,0,0)))
        , cl_helium      (new Atom(2,0.0,Vector3D(0,0,0)))
    {
        BasisSet1::theGlobalCache=new BasisSet1::IntegralsCache_RAM<double>();       
    }
    ~Basis1_RKB_Tests()
    {
        delete cl_hydrogen;
        delete cl_hydrogen_100;
        delete cl_helium;
        delete bs1;
        delete bs2;
        delete BasisSet1::theGlobalCache;
    }
    void Init(nlohmann::json js,BasisSetAtom::Type legacy_type)
    {
        bs1=BasisSet1::Atom::Factory(js,Z);
        bs2=BasisSet1::Atom::Factory(js,Z);
        js["type"]=legacy_type;
        legacy_bs=BasisSetAtom::Factory(js,Z);
    }
    // void InitBSpline6() {Init({{"type",BasisSet1::Atom::Type::BSpline6},{"N", N}, {"rmin", 0.1}, {"rmax", 10}},BasisSetAtom::Type::BSpline6);}
    void InitGaussian() {Init({{"type",BasisSet1::Atom::Type::Gaussian_RKB},{"N", N}, {"emin", 0.1}, {"emax", 10}},BasisSetAtom::Type::Gaussian_RKB);}
    void InitSlater  () {Init({{"type",BasisSet1::Atom::Type::Slater_RKB  },{"N", N}, {"emin", 0.1}, {"emax", 10}},BasisSetAtom::Type::Slater_RKB);}
    using OIBS=BasisSet1::Real_ORKB;
    using Legacy_OIBS=Orbital_RKB_IBS<double>;
    void TestOverlap() const
    {
        auto ibs2=bs2->Iterate<OIBS>().begin();
        auto legacy_ibs=legacy_bs->Iterate<Legacy_OIBS>().begin();
        for (auto ibs1:bs1->Iterate<OIBS>())
        {
            auto& S1=ibs1->Overlap();
            auto& S2=(*ibs2)->Overlap();
            auto& legacyS=(*legacy_ibs)->Overlap();
            EXPECT_EQ(S1,S2);
            EXPECT_EQ(&S1,&S2);
            EXPECT_EQ(S1,legacyS);
            ++ibs2;
            ++legacy_ibs;
        }
    }
    void TestKinetic() const
    {
        auto ibs2=bs2->Iterate<OIBS>().begin();
        auto legacy_ibs=legacy_bs->Iterate<Legacy_OIBS>().begin();
        for (auto ibs1:bs1->Iterate<OIBS>())
        {
            auto& S1=ibs1->Kinetic();
            auto& S2=(*ibs2)->Kinetic();
            auto& legacyS=(*legacy_ibs)->Kinetic();
            EXPECT_EQ(S1,S2);
            EXPECT_EQ(&S1,&S2);
            EXPECT_EQ(S1,legacyS);
            ++ibs2;
            ++legacy_ibs;
        }
    }
    void TestNuclear() const
    {
        auto ibs2=bs2->Iterate<OIBS>().begin();
        auto legacy_ibs=legacy_bs->Iterate<Legacy_OIBS>().begin();
        for (auto ibs1:bs1->Iterate<OIBS>())
        {
            auto& S1=ibs1->Nuclear(cl_hydrogen);
            auto& S2=(*ibs2)->Nuclear(cl_hydrogen);
            auto& S3=(*ibs2)->Nuclear(cl_hydrogen_100);
            auto& S4=(*ibs2)->Nuclear(cl_helium);
            auto& legacyS=(*legacy_ibs)->Nuclear(cl_hydrogen);
            EXPECT_EQ(S1,S2);
            EXPECT_EQ(S1,S3);
            EXPECT_NE(S1,S4);
            EXPECT_EQ(&S1,&S2);
            EXPECT_NE(&S1,&S3);
            EXPECT_NE(&S1,&S4);
            EXPECT_EQ(S1,legacyS);
            ++ibs2;
            ++legacy_ibs;
        }
    }
    void TestRestMass() const
    {
        auto ibs2=bs2->Iterate<OIBS>().begin();
        auto legacy_ibs=legacy_bs->Iterate<Legacy_OIBS>().begin();
        for (auto ibs1:bs1->Iterate<OIBS>())
        {
            auto& S1=ibs1->RestMass();
            auto& S2=(*ibs2)->RestMass();
            auto& legacyS=(*legacy_ibs)->RestMass();
            EXPECT_EQ(S1,S2);
            EXPECT_EQ(&S1,&S2);
            EXPECT_EQ(S1,legacyS);
            ++ibs2;
            ++legacy_ibs;
        }
    }
    
    size_t N,Z;
    Cluster *cl_hydrogen,*cl_hydrogen_100,*cl_helium;
    BasisSet1::Real_BS *bs1,*bs2;
    BasisSet* legacy_bs;
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
